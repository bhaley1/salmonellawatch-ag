"""Tests for src.signals — per-cluster signal derivation."""

import sys
from datetime import date, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import signals
from src.lookups import amr as amr_lookup
from src.lookups import clonal_complex as cc_lookup
from src.lookups import geography as geo_lookup


# ---------------------------------------------------------------------------
# Lookup tables
# ---------------------------------------------------------------------------

def test_geography_canonical_us():
    assert geo_lookup.canonical_country("USA") == "United States"
    assert geo_lookup.canonical_country("usa") == "United States"
    assert geo_lookup.canonical_country("United States") == "United States"
    assert geo_lookup.canonical_country("United States of America") == "United States"


def test_geography_canonical_uk():
    assert geo_lookup.canonical_country("UK") == "United Kingdom"
    assert geo_lookup.canonical_country("Great Britain") == "United Kingdom"


def test_geography_canonical_unknown():
    """Unknown countries pass through unchanged; not coerced."""
    assert geo_lookup.canonical_country("Ruritania") == "Ruritania"
    assert geo_lookup.canonical_country(None) is None
    assert geo_lookup.canonical_country("") is None


def test_geography_continent():
    assert geo_lookup.continent_of("USA") == "North America"
    assert geo_lookup.continent_of("Brazil") == "South America"
    assert geo_lookup.continent_of("Germany") == "Europe"
    assert geo_lookup.continent_of("Japan") == "Asia"
    assert geo_lookup.continent_of("Australia") == "Oceania"
    assert geo_lookup.continent_of("Ruritania") is None


def test_geography_hemisphere():
    assert geo_lookup.hemisphere_of("USA") == "N"
    assert geo_lookup.hemisphere_of("Australia") == "S"
    assert geo_lookup.hemisphere_of("Brazil") == "S"
    assert geo_lookup.hemisphere_of("South Africa") == "S"


def test_is_us_state():
    assert geo_lookup.is_us_state("Maryland")
    assert geo_lookup.is_us_state("Texas")
    assert geo_lookup.is_us_state("DC")
    assert not geo_lookup.is_us_state("Ontario")
    assert not geo_lookup.is_us_state(None)


def test_amr_classes_match():
    classes = amr_lookup.amr_clinical_classes({"tet(M)", "blaTEM-1", "fosX"})
    assert "Tetracyclines" in classes
    # tet(M) only matches tetracyclines; blaTEM-1 is not in our curated list
    assert len(classes) == 1


def test_amr_classes_empty():
    classes = amr_lookup.amr_clinical_classes({"fosX", "lin", "mprF"})
    # None of these are in the first-line resistance list
    assert classes == []


def test_amr_multiple_classes():
    classes = amr_lookup.amr_clinical_classes({"tet(M)", "erm(B)", "dfrG"})
    assert "Tetracyclines" in classes
    assert "Macrolides / lincosamides" in classes
    assert "Trimethoprim / sulfonamides" in classes


def test_st_to_cc_known():
    assert cc_lookup.st_to_cc("ST6") == "CC6"
    assert cc_lookup.st_to_cc("ST1") == "CC1"
    assert cc_lookup.st_to_cc("ST4") == "CC4"
    assert cc_lookup.st_to_cc("ST121") == "CC121"


def test_st_to_cc_unknown():
    assert cc_lookup.st_to_cc("ST9999") is None
    assert cc_lookup.st_to_cc("untypeable") is None
    assert cc_lookup.st_to_cc("novel (closest 6)") is None
    assert cc_lookup.st_to_cc(None) is None


def test_cc_notes():
    notes = cc_lookup.cc_notes("CC4")
    assert notes is not None
    assert "label" in notes
    assert "tooltip" in notes
    assert cc_lookup.cc_notes("CC9999") is None
    assert cc_lookup.cc_notes(None) is None


# ---------------------------------------------------------------------------
# Geographic spread
# ---------------------------------------------------------------------------

def test_spread_single_state():
    members = [
        {"geo_country": "USA", "geo_admin1": "Maryland"},
        {"geo_country": "USA", "geo_admin1": "Maryland"},
    ]
    result = signals.compute_geographic_spread(members)
    assert result["label"] is None
    assert result["n_countries"] == 1
    assert result["n_us_states"] == 1


def test_spread_multi_state():
    members = [
        {"geo_country": "USA", "geo_admin1": "Maryland"},
        {"geo_country": "USA", "geo_admin1": "Texas"},
        {"geo_country": "USA", "geo_admin1": "California"},
    ]
    result = signals.compute_geographic_spread(members)
    assert result["label"] == "multi-state"
    assert result["n_us_states"] == 3


def test_spread_multi_country():
    members = [
        {"geo_country": "USA", "geo_admin1": "Maryland"},
        {"geo_country": "Canada", "geo_admin1": "Ontario"},
        {"geo_country": "Mexico", "geo_admin1": None},
    ]
    result = signals.compute_geographic_spread(members)
    assert result["label"] == "multi-country"
    assert result["n_countries"] == 3


def test_spread_multi_continent():
    members = [
        {"geo_country": "USA", "geo_admin1": "Maryland"},
        {"geo_country": "Germany", "geo_admin1": None},
        {"geo_country": "Japan", "geo_admin1": None},
    ]
    result = signals.compute_geographic_spread(members)
    assert result["label"] == "multi-continent"
    assert result["n_continents"] == 3


def test_spread_multi_hemisphere():
    members = [
        {"geo_country": "USA", "geo_admin1": "Maryland"},
        {"geo_country": "Australia", "geo_admin1": None},
    ]
    result = signals.compute_geographic_spread(members)
    assert result["label"] == "multi-hemisphere"
    assert result["n_hemispheres"] == 2


# ---------------------------------------------------------------------------
# Import signal
# ---------------------------------------------------------------------------

def test_import_no_signal_single_country():
    members = [
        {"source_category": "Human", "geo_country": "USA", "food_origin": None},
        {"source_category": "Food", "geo_country": "USA", "food_origin": None},
    ]
    result = signals.compute_import_signal(members)
    assert result["flag"] is False


def test_import_country_mismatch():
    members = [
        {"source_category": "Human", "geo_country": "USA", "food_origin": None, "pdt_acc": "PDT1"},
        {"source_category": "Food", "geo_country": "Mexico", "food_origin": None, "pdt_acc": "PDT2"},
    ]
    result = signals.compute_import_signal(members)
    assert result["flag"] is True
    assert "country mismatch" in result["evidence"]


def test_import_food_origin_field():
    members = [
        {"source_category": "Food", "geo_country": "USA", "food_origin": "Mexico", "pdt_acc": "PDT1"},
        {"source_category": "Human", "geo_country": "USA", "food_origin": None, "pdt_acc": "PDT2"},
    ]
    result = signals.compute_import_signal(members)
    assert result["flag"] is True
    assert "food_origin field" in result["evidence"]
    assert result["n_food_origin_mismatches"] == 1


# ---------------------------------------------------------------------------
# Travel signal
# ---------------------------------------------------------------------------

def test_travel_signal_detects_keyword():
    members = [
        {"source_category": "Human", "host_disease": "travel-associated listeriosis",
         "isolation_source": None, "pdt_acc": "PDT_T1"},
        {"source_category": "Human", "host_disease": "listeriosis",
         "isolation_source": None, "pdt_acc": "PDT_T2"},
    ]
    result = signals.compute_travel_signal(members)
    assert result["flag"] is True
    assert result["n_annotated"] == 1
    assert "PDT_T1" in result["annotated_pdts"]


def test_travel_signal_negative():
    members = [
        {"source_category": "Human", "host_disease": "listeriosis",
         "isolation_source": "blood", "pdt_acc": "PDT_X"},
    ]
    result = signals.compute_travel_signal(members)
    assert result["flag"] is False
    assert result["n_annotated"] == 0


def test_travel_signal_isolation_source_field():
    """Travel keyword can also live in isolation_source."""
    members = [
        {"source_category": "Human", "host_disease": None,
         "isolation_source": "blood, returning traveler", "pdt_acc": "PDT_T"},
    ]
    result = signals.compute_travel_signal(members)
    assert result["flag"] is True


def test_travel_signal_only_humans():
    """A food isolate with travel-related text shouldn't trigger."""
    members = [
        {"source_category": "Food", "host_disease": "travel-associated",
         "isolation_source": None, "pdt_acc": "PDT_FOOD"},
    ]
    result = signals.compute_travel_signal(members)
    assert result["flag"] is False


def test_travel_signal_case_insensitive():
    members = [
        {"source_category": "Human", "host_disease": "TRAVEL-ASSOCIATED LISTERIOSIS",
         "isolation_source": None, "pdt_acc": "PDT_T"},
    ]
    result = signals.compute_travel_signal(members)
    assert result["flag"] is True


def test_travel_signal_traveller_british():
    members = [
        {"source_category": "Human", "host_disease": "returning traveller",
         "isolation_source": None, "pdt_acc": "PDT_T"},
    ]
    result = signals.compute_travel_signal(members)
    assert result["flag"] is True


# ---------------------------------------------------------------------------
# Submitter diversity
# ---------------------------------------------------------------------------

def test_submitter_diversity_flag_fires_at_five():
    members = [
        {"bioproject_acc": f"PRJNA{i:03d}"} for i in range(5)
    ]
    result = signals.compute_submitter_diversity(members)
    assert result["flag"] is True
    assert result["n_bioprojects"] == 5
    assert result["tier"] == "strong"


def test_submitter_diversity_very_strong_at_ten():
    members = [
        {"bioproject_acc": f"PRJNA{i:03d}"} for i in range(12)
    ]
    result = signals.compute_submitter_diversity(members)
    assert result["flag"] is True
    assert result["tier"] == "very strong"


def test_submitter_diversity_below_threshold():
    members = [
        {"bioproject_acc": "PRJNA001"},
        {"bioproject_acc": "PRJNA001"},  # duplicate — same BioProject
        {"bioproject_acc": "PRJNA002"},
    ]
    result = signals.compute_submitter_diversity(members)
    assert result["flag"] is False
    assert result["n_bioprojects"] == 2


def test_submitter_diversity_handles_missing():
    members = [
        {"bioproject_acc": None},
        {"bioproject_acc": ""},
        {},
    ]
    result = signals.compute_submitter_diversity(members)
    assert result["flag"] is False
    assert result["n_bioprojects"] == 0


# ---------------------------------------------------------------------------
# Human emergence
# ---------------------------------------------------------------------------

def test_human_emergence_classic_pattern():
    """Cluster that was mostly food, now becoming human-dominated."""
    today = date(2026, 5, 21)
    # 15 baseline food isolates from 2-3 years ago, no humans
    baseline = [
        {"source_category": "Food", "collection_date": "2023-06-01"},
    ] * 15
    # 5 recent isolates, 4 of them human
    recent = [
        {"source_category": "Human", "collection_date": "2026-04-15"},
        {"source_category": "Human", "collection_date": "2026-04-10"},
        {"source_category": "Human", "collection_date": "2026-03-15"},
        {"source_category": "Human", "collection_date": "2026-02-01"},
        {"source_category": "Food", "collection_date": "2026-03-20"},
    ]
    result = signals.compute_human_emergence(baseline + recent, today)
    assert result["flag"] is True
    assert result["recent_ratio"] == 0.8
    assert result["baseline_ratio"] == 0.0
    assert result["delta"] == 0.8


def test_human_emergence_stable_cluster_does_not_fire():
    """A cluster that's always been 50/50 doesn't trigger emergence."""
    today = date(2026, 5, 21)
    baseline = []
    for i in range(10):
        baseline.extend([
            {"source_category": "Human", "collection_date": "2023-06-01"},
            {"source_category": "Food", "collection_date": "2023-06-01"},
        ])
    recent = [
        {"source_category": "Human", "collection_date": "2026-04-15"},
        {"source_category": "Food", "collection_date": "2026-04-15"},
        {"source_category": "Human", "collection_date": "2026-03-15"},
    ]
    result = signals.compute_human_emergence(baseline + recent, today)
    assert result["flag"] is False


def test_human_emergence_needs_baseline():
    """Young cluster with no baseline should not fire."""
    today = date(2026, 5, 21)
    members = [
        {"source_category": "Human", "collection_date": "2026-04-15"},
        {"source_category": "Human", "collection_date": "2026-04-10"},
    ]
    result = signals.compute_human_emergence(members, today)
    assert result["flag"] is False
    assert result.get("reason") == "baseline_too_small"


def test_human_emergence_needs_recent_activity():
    """Cluster with strong baseline but no recent activity doesn't fire."""
    today = date(2026, 5, 21)
    members = [
        {"source_category": "Food", "collection_date": "2022-06-01"},
    ] * 15
    result = signals.compute_human_emergence(members, today)
    assert result["flag"] is False
    assert result.get("reason") == "no_recent_activity"


def test_human_emergence_skips_already_human_dominant():
    """If baseline was already mostly human, this isn't 'emergence'."""
    today = date(2026, 5, 21)
    # Baseline already 80% human
    baseline = []
    for i in range(8):
        baseline.append({"source_category": "Human", "collection_date": "2023-06-01"})
    for i in range(2):
        baseline.append({"source_category": "Food", "collection_date": "2023-06-01"})
    # Recent 90% human (small shift, not emergence)
    recent = [
        {"source_category": "Human", "collection_date": "2026-04-15"},
    ] * 9
    recent.append({"source_category": "Food", "collection_date": "2026-04-15"})

    result = signals.compute_human_emergence(baseline + recent, today)
    # Baseline was already >50% human → not "emergence" — should not flag
    assert result["flag"] is False


# ---------------------------------------------------------------------------
# AMR critical
# ---------------------------------------------------------------------------

def test_amr_critical_present():
    amr = {"PDT1": {"tet(M)"}, "PDT2": {"fosX"}}
    result = signals.compute_amr_critical(amr)
    assert result["flag"] is True
    assert "Tetracyclines" in result["classes_present"]


def test_amr_critical_absent():
    """fosX is a typical Listeria intrinsic-resistance gene, not flagged."""
    amr = {"PDT1": {"fosX"}, "PDT2": {"lin", "mprF"}}
    result = signals.compute_amr_critical(amr)
    assert result["flag"] is False
    assert result["classes_present"] == []


def test_amr_critical_multiple_classes():
    amr = {"PDT1": {"tet(M)", "erm(B)", "dfrG"}}
    result = signals.compute_amr_critical(amr)
    assert result["flag"] is True
    assert result["n_classes"] == 3


# ---------------------------------------------------------------------------
# Acceleration
# ---------------------------------------------------------------------------

def test_acceleration_too_young():
    """Clusters younger than baseline_min_days can't accelerate."""
    today = date(2026, 5, 21)
    members = [
        {"target_creation_date": (today - timedelta(days=10)).isoformat()},
        {"target_creation_date": (today - timedelta(days=20)).isoformat()},
    ]
    result = signals.compute_acceleration(members, today)
    assert result["flag"] is False
    assert result.get("reason") == "cluster_too_young"


def test_acceleration_active():
    """Recent rate clearly > baseline → flag."""
    today = date(2026, 5, 21)
    # 200 day cluster: 2 isolates spread over the first 140 days (baseline ~0.014/day),
    # then 5 new isolates in the last 60 days (recent ~0.083/day) → ratio ~6
    baseline_dates = [
        today - timedelta(days=200),
        today - timedelta(days=120),
    ]
    recent_dates = [
        today - timedelta(days=5),
        today - timedelta(days=10),
        today - timedelta(days=20),
        today - timedelta(days=30),
        today - timedelta(days=50),
    ]
    members = [{"target_creation_date": d.isoformat()} for d in baseline_dates + recent_dates]
    result = signals.compute_acceleration(members, today)
    assert result["flag"] is True
    assert result["ratio"] > 2.0


def test_acceleration_stable():
    """Steady rate over a long period → no flag."""
    today = date(2026, 5, 21)
    # 1 isolate every 30 days for 300 days: stable rate
    dates = [today - timedelta(days=i) for i in range(0, 300, 30)]
    members = [{"target_creation_date": d.isoformat()} for d in dates]
    result = signals.compute_acceleration(members, today)
    assert result["flag"] is False


# ---------------------------------------------------------------------------
# IFSAC
# ---------------------------------------------------------------------------

def test_ifsac_top_category():
    members = [
        {"ifsac_category": "Dairy"},
        {"ifsac_category": "Dairy"},
        {"ifsac_category": "Vegetable Row Crops"},
        {"ifsac_category": None},
    ]
    result = signals.compute_ifsac_summary(members)
    assert result["top"]["category"] == "Dairy"
    assert result["top"]["n"] == 2
    assert len(result["categories"]) == 2


def test_ifsac_none():
    members = [
        {"ifsac_category": None},
        {"ifsac_category": ""},
    ]
    result = signals.compute_ifsac_summary(members)
    assert result["top"] is None
    assert result["categories"] == []


# ---------------------------------------------------------------------------
# CC
# ---------------------------------------------------------------------------

def test_cc_known_with_notes():
    result = signals.compute_clonal_complex("ST6")
    assert result["cc"] == "CC6"
    assert result["label"] is not None
    assert result["has_notes"]


def test_cc_known_without_notes():
    """ST that maps to a CC we have no curated notes for."""
    # ST7 maps to CC7 which is in ST_TO_CC but not in CC_NOTES
    result = signals.compute_clonal_complex("ST7")
    assert result["cc"] == "CC7"
    assert result["label"] is None
    assert not result["has_notes"]


def test_cc_unknown_st():
    result = signals.compute_clonal_complex("ST9999")
    assert result["cc"] is None


def test_cc_untypeable():
    result = signals.compute_clonal_complex("untypeable")
    assert result["cc"] is None


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------

def test_compute_all_signals_smoke():
    """Smoke-test that the orchestrator returns the expected shape."""
    members = [
        {"source_category": "Human", "geo_country": "USA", "geo_admin1": "Maryland",
         "food_origin": None, "ifsac_category": None,
         "host_disease": "travel-associated listeriosis", "isolation_source": None,
         "pdt_acc": "PDT_H",
         "target_creation_date": "2026-05-15"},
        {"source_category": "Food", "geo_country": "USA", "geo_admin1": "Texas",
         "food_origin": None, "ifsac_category": "Dairy",
         "host_disease": None, "isolation_source": "milk",
         "pdt_acc": "PDT_F",
         "target_creation_date": "2025-12-01"},
    ]
    result = signals.compute_all_signals(
        members=members,
        amr_genes_per_isolate={"PDT1": {"tet(M)"}},
        mlst_st="ST6",
        today=date(2026, 5, 21),
    )
    assert "geographic_spread" in result
    assert "import_signal" in result
    assert "travel_signal" in result
    assert "submitter_diversity" in result
    assert "human_emergence" in result
    assert "amr_critical" in result
    assert "acceleration" in result
    assert "ifsac" in result
    assert "clonal_complex" in result
    # AMR critical should flag because of tet(M)
    assert result["amr_critical"]["flag"] is True
    # CC should be CC6
    assert result["clonal_complex"]["cc"] == "CC6"
    # Travel signal should flag because of the host_disease annotation
    assert result["travel_signal"]["flag"] is True


# ---------------------------------------------------------------------------
# Test runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import inspect
    mod = sys.modules[__name__]
    tests = [
        (name, fn) for name, fn in inspect.getmembers(mod, inspect.isfunction)
        if name.startswith("test_")
    ]
    passed, failed = 0, 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
            passed += 1
        except Exception as e:
            import traceback
            print(f"  FAIL  {name}: {type(e).__name__}: {e}")
            traceback.print_exc()
            failed += 1
    print()
    print(f"{passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
