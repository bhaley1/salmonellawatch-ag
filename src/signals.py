"""Derive per-cluster signals from cluster member data.

Each function takes the materialized cluster member list (or a small
projection of it) and returns a small dict that gets serialized into
cluster_summary.signals_json and rendered as a chip on the dashboard.

Signals computed here:
  - geographic_spread:    multi-state / multi-country / multi-continent
  - import_signal:        country mismatch between human and nonhuman isolates
  - amr_critical:         first-line resistance gene present
  - acceleration:         recent rate > 2x historical baseline
  - ifsac_summary:        top IFSAC food category in the cluster
  - clonal_complex:       ST → CC mapping, with curated notes
"""

from __future__ import annotations

from collections import Counter
from datetime import date, timedelta

from .lookups import amr as amr_lookup
from .lookups import clonal_complex as cc_lookup
from .lookups import geography as geo_lookup


# ---------------------------------------------------------------------------
# Geographic spread
# ---------------------------------------------------------------------------
def compute_geographic_spread(members: list[dict]) -> dict:
    """Detect cross-state, cross-country, cross-continent, cross-hemisphere spread.

    Returns a dict with the strongest applicable label and supporting counts.
    Only one badge displays per cluster (the strongest signal), but we
    record all four for transparency.
    """
    countries: set[str] = set()
    continents: set[str] = set()
    hemispheres: set[str] = set()
    states: set[str] = set()  # US states specifically (other countries don't expose admin1 as consistently)

    for m in members:
        country = geo_lookup.canonical_country(m.get("geo_country"))
        if country:
            countries.add(country)
            cont = geo_lookup.continent_of(country)
            if cont:
                continents.add(cont)
            hemi = geo_lookup.hemisphere_of(country)
            if hemi:
                hemispheres.add(hemi)
            # Multi-state: only counted for the US (where admin1 = state)
            if country == "United States" and geo_lookup.is_us_state(m.get("geo_admin1")):
                states.add(m["geo_admin1"])

    # Pick the strongest label.
    label = None
    if len(hemispheres) >= 2:
        label = "multi-hemisphere"
    elif len(continents) >= 2:
        label = "multi-continent"
    elif len(countries) >= 2:
        label = "multi-country"
    elif len(states) >= 2:
        label = "multi-state"

    return {
        "label": label,
        "n_countries": len(countries),
        "n_continents": len(continents),
        "n_hemispheres": len(hemispheres),
        "n_us_states": len(states),
    }


# ---------------------------------------------------------------------------
# Import signal — the food vector crossed a border
# ---------------------------------------------------------------------------
def compute_import_signal(members: list[dict]) -> dict:
    """Flag clusters where human cases are in different countries from food isolates.

    Two evidence types:
      (a) Country mismatch: cluster contains human-source isolates in country A
          AND food-source isolates in country B (B != A).
      (b) food_origin field populated and differs from isolation country.

    Returns {flag: bool, evidence: str|None, detail: dict}.
    """
    human_countries: set[str] = set()
    food_countries: set[str] = set()
    food_origin_mismatches: list[dict] = []

    for m in members:
        sc = m.get("source_category")
        country = geo_lookup.canonical_country(m.get("geo_country"))
        food_origin = (m.get("food_origin") or "").strip() or None

        if sc == "Human" and country:
            human_countries.add(country)
        elif sc == "Food" and country:
            food_countries.add(country)

        if food_origin:
            origin_canonical = geo_lookup.canonical_country(food_origin)
            if origin_canonical and country and origin_canonical != country:
                food_origin_mismatches.append({
                    "pdt": m.get("pdt_acc"),
                    "isolation_country": country,
                    "food_origin_country": origin_canonical,
                })

    # Evidence (a): country mismatch between human and food
    country_mismatch = bool(
        human_countries and food_countries and (human_countries - food_countries or food_countries - human_countries)
    )

    flag = country_mismatch or bool(food_origin_mismatches)
    if food_origin_mismatches:
        evidence = "food_origin field"
    elif country_mismatch:
        evidence = "country mismatch between human and food isolates"
    else:
        evidence = None

    return {
        "flag": flag,
        "evidence": evidence,
        "n_food_origin_mismatches": len(food_origin_mismatches),
        "human_countries": sorted(human_countries),
        "food_countries": sorted(food_countries),
    }


# ---------------------------------------------------------------------------
# Travel signal — the patient crossed a border
# ---------------------------------------------------------------------------

# Keywords that, when present in host_disease or isolation_source, suggest
# a travel-associated case. Deliberately conservative: only well-defined
# travel terminology. We don't try to match every possible phrasing because
# false positives are worse than false negatives here — under-detection is
# already a known limitation we document explicitly.
_TRAVEL_KEYWORDS = (
    "travel-associated",
    "travel associated",
    "travel-related",
    "travel related",
    "traveler",
    "traveller",        # British spelling
    "returning traveler",
    "returning traveller",
    "imported case",
    "imported listeriosis",
    "imported infection",
)


def _has_travel_keyword(text: str | None) -> bool:
    if not text:
        return False
    t = text.lower()
    return any(kw in t for kw in _TRAVEL_KEYWORDS)


def compute_travel_signal(members: list[dict]) -> dict:
    """Flag clusters with submitter-annotated travel-associated cases.

    IMPORTANT LIMITATION: NCBI Pathogen Detection has no dedicated travel
    history field. Travel association is only detectable when the submitter
    happened to annotate it in free-text fields (host_disease,
    isolation_source). The vast majority of travel-associated cases are
    NOT annotated this way — they appear as ordinary domestic cases.

    This signal therefore has very high false-negative rate. It is useful
    as a "this cluster has at least one known travel case worth
    investigating" pointer, not as a measure of how travel-associated the
    cluster actually is. The methods page documents this limitation
    explicitly.

    Detection: case-insensitive substring match of curated keywords in
    host_disease or isolation_source for any human cluster member.
    """
    annotated_pdts: list[str] = []
    for m in members:
        # Only human members can have travel-associated illness
        if m.get("source_category") != "Human":
            continue
        hd = m.get("host_disease")
        iso = m.get("isolation_source")
        if _has_travel_keyword(hd) or _has_travel_keyword(iso):
            annotated_pdts.append(m.get("pdt_acc") or "")

    return {
        "flag": len(annotated_pdts) > 0,
        "n_annotated": len(annotated_pdts),
        "annotated_pdts": annotated_pdts[:10],  # cap for storage size
    }


# ---------------------------------------------------------------------------
# Submitter diversity — independent labs detecting the same strain
# ---------------------------------------------------------------------------
def compute_submitter_diversity(members: list[dict]) -> dict:
    """Count distinct BioProjects in the cluster.

    BioProject (PRJNA...) is NCBI's grouping for "a set of related research
    or surveillance data submitted by one organization." It's an imperfect
    proxy for "independent labs detecting this strain" — large agencies use
    one BioProject for years of submissions, and a small lab might use
    multiple BioProjects for one study — but in aggregate it's a useful
    submitter-diversity signal.

    When many independent BioProjects detect genomically-related isolates,
    that's much harder to explain away as cohort effects or one-lab bias
    than when the same BioProject keeps re-submitting closely-related
    isolates. We flag clusters with ≥5 distinct BioProjects as a strong
    submitter-diversity signal.
    """
    bioprojects: set[str] = set()
    for m in members:
        bp = (m.get("bioproject_acc") or "").strip()
        if bp:
            bioprojects.add(bp)

    n = len(bioprojects)
    # Tier thresholds:
    #   ≥10 BioProjects: very strong signal (label)
    #   ≥5 BioProjects:  strong signal (flag fires)
    #   2-4 BioProjects: minor signal (counted but no chip)
    #   1 or 0:          single-source cluster (no chip)
    flag = n >= 5
    tier = None
    if n >= 10:
        tier = "very strong"
    elif n >= 5:
        tier = "strong"
    return {
        "flag": flag,
        "n_bioprojects": n,
        "tier": tier,
    }


# ---------------------------------------------------------------------------
# Human-emergence signal — cluster transitioning from food-dominant to
# human-dominant is the canonical "outbreak emerging" pattern
# ---------------------------------------------------------------------------
def compute_human_emergence(
    members: list[dict],
    today: date,
    recent_window_days: int = 365,
) -> dict:
    """Detect clusters where human cases are appearing in a previously
    food/environment-dominated cluster.

    The epidemiological logic: many Listeria clusters start as food or
    environmental contamination signals (a processing facility persistently
    sheds the strain into food production for years). When human cases
    start to appear in such a cluster, that's the signal that the
    contamination has finally caused detected disease — the canonical
    "outbreak emerging" pattern that surveillance professionals watch for.

    Detection:
      - Compute the human:total ratio over the recent window (last 365 days)
        and over the cluster's prior baseline (all members older than that).
      - Flag when recent ratio is meaningfully higher than baseline AND the
        cluster has ≥10 prior isolates (so we have a real baseline).
      - "Meaningfully higher" means: recent_ratio ≥ baseline_ratio + 0.20
        AND recent_ratio ≥ 0.30 (so a cluster that goes from 5% to 25%
        human triggers, but a cluster that goes from 5% to 6% does not).

    Returns {flag, recent_ratio, baseline_ratio, n_recent_human, n_recent_total,
             reason}.
    """
    cutoff = (today - timedelta(days=recent_window_days)).isoformat()

    n_recent_total = 0
    n_recent_human = 0
    n_baseline_total = 0
    n_baseline_human = 0

    for m in members:
        cd = m.get("collection_date")
        if not cd:
            continue
        is_recent = str(cd)[:10] >= cutoff
        is_human = m.get("source_category") == "Human"
        if is_recent:
            n_recent_total += 1
            if is_human:
                n_recent_human += 1
        else:
            n_baseline_total += 1
            if is_human:
                n_baseline_human += 1

    # Need a real baseline
    if n_baseline_total < 10:
        return {
            "flag": False,
            "n_recent_human": n_recent_human,
            "n_recent_total": n_recent_total,
            "n_baseline_total": n_baseline_total,
            "reason": "baseline_too_small",
        }
    # Need recent activity
    if n_recent_total < 3:
        return {
            "flag": False,
            "n_recent_human": n_recent_human,
            "n_recent_total": n_recent_total,
            "n_baseline_total": n_baseline_total,
            "reason": "no_recent_activity",
        }

    recent_ratio = n_recent_human / n_recent_total
    baseline_ratio = n_baseline_human / n_baseline_total

    # Flag conditions
    delta = recent_ratio - baseline_ratio
    flag = (delta >= 0.20) and (recent_ratio >= 0.30) and (baseline_ratio < 0.50)

    return {
        "flag": flag,
        "recent_ratio": round(recent_ratio, 3),
        "baseline_ratio": round(baseline_ratio, 3),
        "delta": round(delta, 3),
        "n_recent_human": n_recent_human,
        "n_recent_total": n_recent_total,
        "n_baseline_total": n_baseline_total,
    }


# ---------------------------------------------------------------------------
# AMR clinical importance
# ---------------------------------------------------------------------------
def compute_amr_critical(amr_genes_per_isolate: dict[str, set[str]]) -> dict:
    """Flag clusters carrying first-line resistance genes.

    Input: {pdt_acc: {gene_symbol, gene_symbol, ...}} for all cluster members.
    Returns {flag, classes_present, gene_count, ...}.
    """
    all_genes: set[str] = set()
    for gset in amr_genes_per_isolate.values():
        all_genes.update(gset)
    classes = amr_lookup.amr_clinical_classes(all_genes)

    # Per-class details: which genes and how many isolates
    class_details: dict[str, dict] = {}
    for cls in classes:
        class_genes_set = amr_lookup.FIRST_LINE_RESISTANCE_GENES[cls]
        hits: dict[str, int] = {}
        for pdt, isolate_genes in amr_genes_per_isolate.items():
            present = isolate_genes & class_genes_set
            for g in present:
                hits[g] = hits.get(g, 0) + 1
        if hits:
            class_details[cls] = hits

    return {
        "flag": len(classes) > 0,
        "classes_present": classes,
        "class_details": class_details,
        "n_classes": len(classes),
    }


# ---------------------------------------------------------------------------
# Acceleration
# ---------------------------------------------------------------------------
def compute_acceleration(
    members: list[dict],
    today: date,
    recent_window_days: int = 60,
    baseline_min_days: int = 180,
) -> dict:
    """Detect whether the cluster is accelerating.

    Compare per-day addition rate over the recent window to the per-day rate
    over the cluster's full history. If the cluster history is shorter than
    baseline_min_days, we can't form a meaningful baseline.

    Returns {flag: bool, recent_rate: float, baseline_rate: float, ratio: float|None}.
    """
    # Use target_creation_date — when the isolate became visible to NCBI's
    # surveillance. That's the right denominator for "is this cluster
    # picking up speed in the surveillance signal."
    dates: list[date] = []
    for m in members:
        td = m.get("target_creation_date")
        if not td:
            continue
        if isinstance(td, str):
            try:
                td = date.fromisoformat(td[:10])
            except ValueError:
                continue
        dates.append(td)

    if not dates:
        return {"flag": False, "recent_rate": 0.0, "baseline_rate": 0.0, "ratio": None}

    earliest = min(dates)
    cluster_age_days = (today - earliest).days

    if cluster_age_days < baseline_min_days:
        # Cluster is too young to have a meaningful baseline
        return {
            "flag": False,
            "recent_rate": 0.0,
            "baseline_rate": 0.0,
            "ratio": None,
            "reason": "cluster_too_young",
        }

    cutoff = today - timedelta(days=recent_window_days)
    n_recent = sum(1 for d in dates if d >= cutoff)
    n_baseline = len(dates) - n_recent
    baseline_days = cluster_age_days - recent_window_days

    recent_rate = n_recent / recent_window_days if recent_window_days > 0 else 0.0
    baseline_rate = n_baseline / baseline_days if baseline_days > 0 else 0.0

    if baseline_rate <= 0:
        # No prior history, or all isolates are in the recent window
        # — treat that as "newly active" rather than "accelerating"
        return {
            "flag": False,
            "recent_rate": recent_rate,
            "baseline_rate": 0.0,
            "ratio": None,
            "reason": "no_prior_baseline",
        }

    ratio = recent_rate / baseline_rate
    flag = ratio >= 2.0 and n_recent >= 2  # at least 2 recent isolates for the signal to be meaningful

    return {
        "flag": flag,
        "recent_rate": round(recent_rate, 4),
        "baseline_rate": round(baseline_rate, 4),
        "ratio": round(ratio, 2),
        "n_recent": n_recent,
        "n_baseline": n_baseline,
    }


# ---------------------------------------------------------------------------
# IFSAC food categorization
# ---------------------------------------------------------------------------
def compute_ifsac_summary(members: list[dict]) -> dict:
    """Summarize IFSAC categories present in the cluster.

    IFSAC = Interagency Food Safety Analytics Collaboration. Categories are
    standardized food groupings (e.g., "Dairy", "Vegetable Row Crops",
    "Poultry"). NCBI's IFSAC_category column when populated is the most
    authoritative food categorization available.

    Counts NONHUMAN isolates only. NCBI sometimes records IFSAC categories
    for clinical isolates ("the patient ate X") but that's a different
    epidemiological semantic from "isolate sourced from food X" and would
    inflate the food signal.
    """
    counter: Counter[str] = Counter()
    for m in members:
        # Skip human isolates — IFSAC on humans is about reported food
        # exposure, not isolate source
        if m.get("source_category") == "Human":
            continue
        ifsac = m.get("ifsac_category")
        if ifsac and ifsac.strip() and ifsac.upper() != "NULL":
            counter[ifsac.strip()] += 1

    if not counter:
        return {"top": None, "categories": []}

    top = counter.most_common(1)[0]
    return {
        "top": {"category": top[0], "n": top[1]},
        "categories": [{"category": c, "n": n} for c, n in counter.most_common(5)],
        "n_total_ifsac": sum(counter.values()),
    }


# ---------------------------------------------------------------------------
# Clonal complex
# ---------------------------------------------------------------------------
def compute_clonal_complex(mlst_st: str | None) -> dict:
    """Derive CC from MLST ST and attach curated notes if any."""
    cc = cc_lookup.st_to_cc(mlst_st)
    notes = cc_lookup.cc_notes(cc) if cc else None
    return {
        "cc": cc,
        "label": notes["label"] if notes else None,
        "tooltip": notes["tooltip"] if notes else None,
        "has_notes": notes is not None,
    }


# ---------------------------------------------------------------------------
# Orchestrator: compute everything for one cluster
# ---------------------------------------------------------------------------
def compute_all_signals(
    members: list[dict],
    amr_genes_per_isolate: dict[str, set[str]],
    mlst_st: str | None,
    today: date,
    recent_window_days: int = 60,
) -> dict:
    """Compute every signal for a cluster, returning a single dict."""
    return {
        "geographic_spread": compute_geographic_spread(members),
        "import_signal": compute_import_signal(members),
        "travel_signal": compute_travel_signal(members),
        "submitter_diversity": compute_submitter_diversity(members),
        "human_emergence": compute_human_emergence(members, today),
        "amr_critical": compute_amr_critical(amr_genes_per_isolate),
        "acceleration": compute_acceleration(members, today, recent_window_days),
        "ifsac": compute_ifsac_summary(members),
        "clonal_complex": compute_clonal_complex(mlst_st),
    }
