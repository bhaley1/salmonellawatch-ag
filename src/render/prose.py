"""Prose helpers — turn raw numbers into human-readable summaries.

Kept small and pure. No DB access, no I/O.
"""

from __future__ import annotations


def format_temporal_span(days: int | None) -> str:
    """Human-readable temporal span with epidemiological framing.

    Returns labels that hint at what the span *means* for Listeria
    surveillance interpretation:
      - Very short (days) → "isolated event"
      - Weeks → "recent activity"
      - Months → "extended outbreak window"
      - Years → "persistent / recurring source"
    """
    if days is None:
        return "—"
    if days < 1:
        return "single day"
    if days < 14:
        return f"{days} days"
    if days < 60:
        weeks = round(days / 7)
        return f"{weeks} weeks"
    if days < 365 * 2:
        months = round(days / 30.4)
        return f"{months} months"
    years = round(days / 365.25, 1)
    if years == int(years):
        return f"{int(years)} years"
    return f"{years} years"


def span_interpretation(days: int | None) -> str | None:
    """Short interpretive phrase for the span. None if not interpretable.

    These framings reflect Listeria's specific epidemiology: short spans
    suggest discrete outbreaks; long spans suggest persistent contamination
    (especially in food processing environments).
    """
    if days is None:
        return None
    if days < 60:
        return "recent activity"
    if days < 365:
        return "extended outbreak window"
    if days < 365 * 3:
        return "multi-year cluster"
    return "long-persistent (possible environmental reservoir)"


def format_geography(admin1_data: dict, max_countries: int = 6, max_admin1_per_country: int = 4) -> str:
    """Format admin1 breakdown into a compact human-readable string.

    Three problems this fixes:
      1. NCBI sometimes records the country name in both the country and
         admin1 fields (e.g. "United Kingdom: United Kingdom"). When
         admin1 == country, treat as unspecified instead of as a redundant
         admin1 entry.
      2. We previously listed every country's "N unspecified" separately,
         which produced very long lines. Now we just suffix per country.
      3. Long country lists are capped at max_countries with "+N more"
         indicating the tail.

    Input shape: {"by_country": {country: [{admin1, n}, ...]}, "unspecified": {country: n}}
    """
    if not admin1_data:
        return "—"

    # Defensive copies so we can mutate
    by_country: dict = dict(admin1_data.get("by_country") or {})
    unspecified: dict[str, int] = dict(admin1_data.get("unspecified") or {})

    # Fix (1): an admin1 that equals its country name is really unspecified.
    # Move that count from by_country into unspecified.
    cleaned_by_country: dict[str, list[dict]] = {}
    for country, admin1_list in by_country.items():
        kept_admin1: list[dict] = []
        for entry in admin1_list:
            a = (entry.get("admin1") or "").strip()
            if not a or a.lower() == country.lower():
                # Redundant — fold into unspecified
                unspecified[country] = unspecified.get(country, 0) + int(entry.get("n", 0))
            else:
                kept_admin1.append(entry)
        if kept_admin1:
            cleaned_by_country[country] = kept_admin1

    # Build per-country totals so we can rank
    country_totals: dict[str, int] = {}
    for country, admin1_list in cleaned_by_country.items():
        country_totals[country] = sum(int(a.get("n", 0)) for a in admin1_list)
    for country, n in unspecified.items():
        country_totals[country] = country_totals.get(country, 0) + n

    if not country_totals:
        return "—"

    sorted_countries = sorted(country_totals.items(), key=lambda x: -x[1])
    visible = sorted_countries[:max_countries]
    hidden_count = sum(n for _, n in sorted_countries[max_countries:])
    hidden_countries = len(sorted_countries) - len(visible)

    parts: list[str] = []
    for country, total in visible:
        admin1_list = cleaned_by_country.get(country, [])[:max_admin1_per_country]
        admin1_str = ", ".join(f"{a['admin1']} ({a['n']})" for a in admin1_list)
        remainder_admin1 = len(cleaned_by_country.get(country, [])) - max_admin1_per_country
        unspec_n = unspecified.get(country, 0)

        if admin1_str:
            parts_inner = [admin1_str]
            if remainder_admin1 > 0:
                parts_inner.append(f"+{remainder_admin1} more")
            if unspec_n:
                parts_inner.append(f"{unspec_n} unspecified")
            parts.append(f"{country} ({total}): {', '.join(parts_inner)}")
        else:
            # No admin1, all unspecified
            parts.append(f"{country} ({total}, location not specified)")

    out = " · ".join(parts)
    if hidden_countries > 0:
        out += f" · +{hidden_countries} more countries ({hidden_count} isolates)"

    return out


def format_oldest_isolate(d: dict | None) -> str:
    """Format an oldest-isolate dict for inline display.

    Input dict shape (from cluster_summary.oldest_*_json):
      {pdt, biosample, date, date_raw, geo, geo_country, geo_admin1, source, source_category}
    """
    if not d:
        return "—"

    date_str = d.get("date_raw") or d.get("date") or "date unknown"
    pdt = d.get("pdt") or "?"
    geo = d.get("geo") or "location unknown"
    source = d.get("source") or "source unknown"
    return f"{date_str} · {geo} · {source} ({pdt})"


def ncbi_link_for(pdt: str | None, biosample: str | None = None) -> str:
    """Pick the best NCBI URL for a PDT accession.

    Prefer the BioSample page (richer metadata) when we have the BioSample
    accession; fall back to NCBI Pathogen Detection's isolate browser
    fragment URL otherwise.
    """
    if biosample:
        return f"https://www.ncbi.nlm.nih.gov/biosample/{biosample}"
    if pdt:
        return f"https://www.ncbi.nlm.nih.gov/pathogens/isolates/#{pdt}"
    return "#"