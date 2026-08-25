"""End-to-end pipeline test using synthetic NCBI fixtures.

Exercises: parse → upsert → summarize → query. No network calls.
Each test runs against a fresh temp SQLite to ensure isolation.
"""

import json
import sys
import tempfile
from datetime import date, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import config, db
from src.ingest import parse, summarize, upsert
from src.render import queries


FIXTURES = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# Test infrastructure
# ---------------------------------------------------------------------------

def fresh_db():
    """Return a temp Path for a new SQLite DB, initialized with schema."""
    fd = tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False)
    fd.close()
    db_path = Path(fd.name)
    db_path.unlink()  # delete; init_db will recreate
    db.init_db(db_path)
    return db_path


def load_fixture_into_db(db_path: Path, pathogen: str = "Listeria"):
    """Run the full parse+upsert+summarize chain against fixture files."""
    pds_map = parse.load_pds_map(FIXTURES / "listeria_cluster_list.tsv")
    isolates = list(parse.iter_isolates(
        FIXTURES / "listeria_metadata.tsv",
        pathogen,
        pds_map,
        pdg_release="PDG_TEST.1",
    ))
    amr_rows = list(parse.iter_amr_rows(FIXTURES / "listeria_amr.tsv"))

    with db.connect(db_path) as conn:
        upsert.upsert_isolates(conn, isolates)
        upsert.upsert_amr(conn, amr_rows)
        upsert.record_release(
            conn, pathogen=pathogen, pdg_release="PDG_TEST.1",
            metadata_url="x", metadata_bytes=1,
            cluster_list_url="x", cluster_list_bytes=1,
            amr_url="x", amr_bytes=1,
        )
        # "today" of 2026-05-21 with a 60-day window catches isolates added
        # from 2026-03-22 onward. Fixture human dates: 5/15, 5/12, 5/8 (inside),
        # 12/25/2025 (outside), and 5/19 (inside).
        summarize.materialize_cluster_summary(
            conn,
            today=date(2026, 5, 21),
        )

    return isolates, amr_rows


# ---------------------------------------------------------------------------
# Parse-layer tests
# ---------------------------------------------------------------------------

def test_parse_pds_map():
    pm = parse.load_pds_map(FIXTURES / "listeria_cluster_list.tsv")
    assert len(pm) == 9
    assert pm["PDT_HUMAN_001.1"] == "PDS_LIST_001.1"
    assert pm["PDT_OTHER_001.1"] == "PDS_LIST_002.1"
    assert "PDT_NOCLUSTER.1" not in pm


def test_parse_isolates_count():
    pm = parse.load_pds_map(FIXTURES / "listeria_cluster_list.tsv")
    isos = list(parse.iter_isolates(
        FIXTURES / "listeria_metadata.tsv",
        "Listeria",
        pm,
        pdg_release="X.1",
    ))
    assert len(isos) == 10  # All 10 rows parsed (including non-clustered)


def test_parse_isolates_source_category():
    """Verify source_category inference is correct for the synthetic fixture."""
    pm = parse.load_pds_map(FIXTURES / "listeria_cluster_list.tsv")
    isos = {i.pdt_acc: i for i in parse.iter_isolates(
        FIXTURES / "listeria_metadata.tsv", "Listeria", pm, "X.1",
    )}
    assert isos["PDT_HUMAN_001.1"].source_category == "Human"
    assert isos["PDT_HUMAN_002.1"].source_category == "Human"
    assert isos["PDT_FOOD_001.1"].source_category == "Food"
    assert isos["PDT_FOOD_002.1"].source_category == "Food"
    # food processing facility drain — Environment wins over Food
    assert isos["PDT_ENV_001.1"].source_category == "Environment"
    # Cattle/milk — Food wins (milk is a food signal in surveillance)
    assert isos["PDT_OTHER_002.1"].source_category == "Food"
    assert isos["PDT_NOCLUSTER.1"].source_category == "Human"


def test_parse_isolates_pds_assignment():
    """Verify cluster assignments propagate from cluster_list."""
    pm = parse.load_pds_map(FIXTURES / "listeria_cluster_list.tsv")
    isos = {i.pdt_acc: i for i in parse.iter_isolates(
        FIXTURES / "listeria_metadata.tsv", "Listeria", pm, "X.1",
    )}
    assert isos["PDT_HUMAN_001.1"].pds_acc == "PDS_LIST_001.1"
    assert isos["PDT_OTHER_001.1"].pds_acc == "PDS_LIST_002.1"
    assert isos["PDT_NOCLUSTER.1"].pds_acc is None


def test_parse_geo_parsing():
    pm = parse.load_pds_map(FIXTURES / "listeria_cluster_list.tsv")
    isos = {i.pdt_acc: i for i in parse.iter_isolates(
        FIXTURES / "listeria_metadata.tsv", "Listeria", pm, "X.1",
    )}
    a = isos["PDT_HUMAN_001.1"]
    assert a.geo_country == "USA"
    assert a.geo_admin1 == "Maryland"


def test_parse_dates():
    """Both collection and target_creation should parse, both raw + parsed."""
    pm = parse.load_pds_map(FIXTURES / "listeria_cluster_list.tsv")
    isos = {i.pdt_acc: i for i in parse.iter_isolates(
        FIXTURES / "listeria_metadata.tsv", "Listeria", pm, "X.1",
    )}
    a = isos["PDT_HUMAN_001.1"]
    assert a.collection_date == date(2026, 4, 15)
    assert a.collection_date_raw == "2026-04-15"
    assert a.target_creation_date == date(2026, 5, 15)


# ---------------------------------------------------------------------------
# Upsert-layer tests
# ---------------------------------------------------------------------------

def test_upsert_idempotent():
    """Running ingest twice should not duplicate rows."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        n1 = conn.execute("SELECT COUNT(*) FROM isolates").fetchone()[0]
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        n2 = conn.execute("SELECT COUNT(*) FROM isolates").fetchone()[0]
    assert n1 == n2 == 10


def test_upsert_amr():
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        n = conn.execute("SELECT COUNT(*) FROM isolate_amr").fetchone()[0]
        assert n == 11
        # PDT_HUMAN_001.1 should have 3 gene calls
        n_human1 = conn.execute(
            "SELECT COUNT(*) FROM isolate_amr WHERE pdt_acc = ?",
            ("PDT_HUMAN_001.1",),
        ).fetchone()[0]
        assert n_human1 == 3


# ---------------------------------------------------------------------------
# Summarize-layer tests
# ---------------------------------------------------------------------------

def test_summarize_basic_counts():
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        # Two clusters in the fixture: PDS_LIST_001.1 (7 members) and PDS_LIST_002.1 (2 members)
        rows = conn.execute(
            "SELECT pds_acc, n_total, n_human, n_nonhuman, n_food, n_animal, n_environment "
            "FROM cluster_summary ORDER BY pds_acc"
        ).fetchall()
    assert len(rows) == 2
    c1 = dict(rows[0])
    assert c1["pds_acc"] == "PDS_LIST_001.1"
    assert c1["n_total"] == 7
    assert c1["n_human"] == 4
    assert c1["n_nonhuman"] == 3
    assert c1["n_food"] == 2
    assert c1["n_environment"] == 1

    c2 = dict(rows[1])
    assert c2["pds_acc"] == "PDS_LIST_002.1"
    assert c2["n_total"] == 2
    assert c2["n_human"] == 1
    # PDT_OTHER_002.1 is Cattle/milk → Food (milk is a food signal)
    assert c2["n_food"] == 1
    assert c2["n_animal"] == 0


def test_summarize_recent_window():
    """Recent windows now filter on collection_date, not target_creation_date.

    With today = 2026-05-21:
      cutoff_60 = 2026-03-22 — captures isolates collected from late March on
      cutoff_30 = 2026-04-21 — captures isolates collected from late April on
      cutoff_15 = 2026-05-06 — captures isolates collected from early May on

    Fixture human collection dates:
      H001 collected 2026-04-15 → in 60d only
      H002 collected 2026-04-10 → in 60d only
      H003 collected 2026-03-28 → in 60d only (March 28 > March 22)
      H004 collected 2025-12-10 → outside all windows
      OTHER_001 collected 2026-05-10 → in 60d, 30d, AND 15d
    """
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT pds_acc, new_humans_in_window, new_humans_30d, new_humans_15d "
            "FROM cluster_summary ORDER BY pds_acc"
        ).fetchall()
    by_pds = {r["pds_acc"]: dict(r) for r in rows}
    # 60d counts: cluster 001 has 3 humans in window, cluster 002 has 1
    assert by_pds["PDS_LIST_001.1"]["new_humans_in_window"] == 3
    assert by_pds["PDS_LIST_002.1"]["new_humans_in_window"] == 1
    # 30d counts: none of cluster 001's humans, but cluster 002's H_OTHER_001 (5/10)
    assert by_pds["PDS_LIST_001.1"]["new_humans_30d"] == 0
    assert by_pds["PDS_LIST_002.1"]["new_humans_30d"] == 1
    # 15d counts: same, only OTHER_001 (5/10 > 5/06)
    assert by_pds["PDS_LIST_001.1"]["new_humans_15d"] == 0
    assert by_pds["PDS_LIST_002.1"]["new_humans_15d"] == 1


def test_summarize_recent_excludes_backlogged_deposit():
    """A human collected long ago but newly deposited should NOT count as a
    new case. The filter must use collection_date, not target_creation_date.
    """
    # H004 in the fixture was collected 2025-12-10 but added 2025-12-25.
    # Both dates are outside the windows from today=2026-05-21, so this
    # test is implicit in the count above. But verify by checking that
    # H004 doesn't appear in the new_humans payload of either cluster.
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT new_humans_in_window_pdts_json FROM cluster_summary"
        ).fetchall()
    for r in rows:
        pdts = json.loads(r["new_humans_in_window_pdts_json"] or "[]")
        assert "PDT_HUMAN_004.1" not in pdts


def test_summarize_countries_json():
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        r = conn.execute(
            "SELECT countries_json FROM cluster_summary WHERE pds_acc = ?",
            ("PDS_LIST_001.1",),
        ).fetchone()
    countries = json.loads(r["countries_json"])
    # All 7 PDS_LIST_001 members are in USA
    assert countries[0]["country"] == "USA"
    assert countries[0]["n"] == 7


def test_summarize_source_summary():
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        r = conn.execute(
            "SELECT source_summary_json FROM cluster_summary WHERE pds_acc = ?",
            ("PDS_LIST_001.1",),
        ).fetchone()
    sources = json.loads(r["source_summary_json"])
    # Should have entries for Food, Environment, and Human categories
    # (source summary now includes all categories, not just nonhuman)
    cats = {s["category"] for s in sources}
    assert "Food" in cats
    assert "Environment" in cats
    assert "Human" in cats  # humans now included


# ---------------------------------------------------------------------------
# Render-layer tests (against queries module)
# ---------------------------------------------------------------------------

def test_query_recent_activity():
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        active = queries.get_recent_activity_clusters(conn)
    # Both clusters have recent activity
    assert len(active) == 2
    # Sorted by new_humans_in_window DESC, so cluster 001 (3 new) comes first
    assert active[0]["pds_acc"] == "PDS_LIST_001.1"
    assert active[0]["new_humans_in_window"] == 3
    assert active[1]["pds_acc"] == "PDS_LIST_002.1"
    assert active[1]["new_humans_in_window"] == 1


def test_query_totals():
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        t = queries.get_totals(conn)
    assert t["n_clusters"] == 2
    assert t["n_human_clusters"] == 2
    assert t["n_mixed"] == 2
    assert t["n_active"] == 2
    assert t["n_new_humans_window"] == 4


def test_summarize_temporal_span():
    """The fixture's earliest collection date is 2025-08-20 (Cattle/milk in
    cluster 002), and latest is 2026-05-10 (PDT_OTHER_001.1). Cluster 001
    spans 2025-10-15 (env) → 2026-05-10 (no wait, 2026-04-15 for H001).
    Verify the temporal_span_days makes sense."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT pds_acc, earliest_collection_date, latest_collection_date, temporal_span_days "
            "FROM cluster_summary ORDER BY pds_acc"
        ).fetchall()
    by_pds = {r["pds_acc"]: dict(r) for r in rows}

    c1 = by_pds["PDS_LIST_001.1"]
    # Cluster 001 dates: 10/15/25 (env), 11/1/25 (food x2), 12/10/25, 3/28/26, 4/10/26, 4/15/26
    # earliest = 2025-10-15, latest = 2026-04-15 → ~182 days
    assert c1["earliest_collection_date"] == "2025-10-15"
    assert c1["latest_collection_date"] == "2026-04-15"
    assert c1["temporal_span_days"] is not None
    assert 180 <= c1["temporal_span_days"] <= 185

    c2 = by_pds["PDS_LIST_002.1"]
    # Cluster 002 dates: 2025-08-20 (Cattle), 2026-05-10 (clinical)
    # span ~263 days
    assert c2["earliest_collection_date"] == "2025-08-20"
    assert c2["latest_collection_date"] == "2026-05-10"
    assert 260 <= c2["temporal_span_days"] <= 265


def test_prose_temporal_span_labels():
    """Verify the human-readable labels for temporal spans."""
    from src.render import prose
    assert prose.format_temporal_span(None) == "—"
    assert prose.format_temporal_span(0) == "single day"
    assert prose.format_temporal_span(5) == "5 days"
    assert prose.format_temporal_span(14) == "2 weeks"
    assert prose.format_temporal_span(45) == "6 weeks"
    assert prose.format_temporal_span(90) == "3 months"
    assert prose.format_temporal_span(365) == "12 months"
    assert prose.format_temporal_span(900) == "2.5 years"
    assert prose.format_temporal_span(365 * 5) == "5 years"


def test_prose_span_interpretation():
    from src.render import prose
    assert prose.span_interpretation(None) is None
    assert prose.span_interpretation(30) == "recent activity"
    assert prose.span_interpretation(180) == "extended outbreak window"
    assert prose.span_interpretation(730) == "multi-year cluster"
    assert prose.span_interpretation(365 * 5) == "long-persistent (possible environmental reservoir)"


def test_query_pathogen_counts():
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        pcs = queries.get_pathogen_counts(conn)
    assert len(pcs) == 1
    assert pcs[0]["pathogen"] == "Listeria"
    assert pcs[0]["n_clusters"] == 2
    assert pcs[0]["n_active"] == 2


def test_recent_humans_details():
    """The new_humans_dates JSON should carry per-isolate detail."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        active = queries.get_recent_activity_clusters(conn)
    c1 = next(c for c in active if c["pds_acc"] == "PDS_LIST_001.1")
    assert len(c1["new_humans_dates"]) == 3
    # Sorted by most-recent collection_date first
    first = c1["new_humans_dates"][0]
    assert first["pdt"] == "PDT_HUMAN_001.1"  # collected 2026-04-15 — most recent
    assert first["geo"] == "USA: Maryland"
    # date_added is still in the payload for transparency, but template won't show it
    assert first["date_added"] == "2026-05-15"
    # biosample is now included for the NCBI link
    assert first["biosample"] == "SAMN12345001"


def test_summarize_oldest_isolates():
    """The cluster summary should include three oldest-isolate signatures."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        r = conn.execute(
            "SELECT oldest_isolate_json, oldest_human_json, oldest_nonhuman_json "
            "FROM cluster_summary WHERE pds_acc = ?",
            ("PDS_LIST_001.1",),
        ).fetchone()

    oldest_any = json.loads(r["oldest_isolate_json"])
    oldest_human = json.loads(r["oldest_human_json"])
    oldest_nonhuman = json.loads(r["oldest_nonhuman_json"])

    # In PDS_LIST_001.1, the oldest isolate of any kind is the env sample
    # PDT_ENV_001.1 collected 2025-10-15
    assert oldest_any is not None
    assert oldest_any["pdt"] == "PDT_ENV_001.1"
    assert oldest_any["date"] == "2025-10-15"
    # biosample is in the pack for NCBI linking
    assert "biosample" in oldest_any

    # Oldest human in PDS_LIST_001.1 is PDT_HUMAN_004.1 collected 2025-12-10
    assert oldest_human is not None
    assert oldest_human["pdt"] == "PDT_HUMAN_004.1"
    assert oldest_human["date"] == "2025-12-10"

    # Oldest nonhuman is the env sample again
    assert oldest_nonhuman is not None
    assert oldest_nonhuman["pdt"] == "PDT_ENV_001.1"


def test_summarize_deposit_lag():
    """The summarizer should compute median + mean deposit lag per cluster."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT pds_acc, deposit_lag_median, deposit_lag_mean "
            "FROM cluster_summary ORDER BY pds_acc"
        ).fetchall()
    by_pds = {r["pds_acc"]: dict(r) for r in rows}
    # PDS_LIST_001.1 isolates have lags around 20-30 days
    # (e.g. PDT_FOOD_001: collected 2025-11-01, deposited 2025-11-20 = 19 days)
    assert by_pds["PDS_LIST_001.1"]["deposit_lag_median"] is not None
    assert by_pds["PDS_LIST_001.1"]["deposit_lag_median"] > 0
    assert by_pds["PDS_LIST_001.1"]["deposit_lag_median"] < 60


def test_summarize_host_summary():
    """The summarizer should aggregate animal hosts."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT pds_acc, host_summary_json "
            "FROM cluster_summary ORDER BY pds_acc"
        ).fetchall()
    by_pds = {r["pds_acc"]: r for r in rows}
    # PDS_LIST_002.1 has PDT_OTHER_002 with host=Bos taurus, no other animals
    hs_002 = json.loads(by_pds["PDS_LIST_002.1"]["host_summary_json"] or "[]")
    hosts_002 = {h["host"]: h["n"] for h in hs_002}
    assert "bos taurus" in hosts_002


def test_summarize_bioproject_diversity():
    """The submitter_diversity signal should count distinct BioProjects."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT pds_acc, signals_json FROM cluster_summary ORDER BY pds_acc"
        ).fetchall()
    by_pds = {r["pds_acc"]: r for r in rows}
    sigs_001 = json.loads(by_pds["PDS_LIST_001.1"]["signals_json"])
    # Cluster 001 has 6 distinct BioProjects (PRJNA001, PRJNA002, PRJNA003,
    # PRJNA010, PRJNA011, PRJNA012) → flag should fire
    assert sigs_001["submitter_diversity"]["n_bioprojects"] >= 5
    assert sigs_001["submitter_diversity"]["flag"] is True


def test_summarize_latest_assembly():
    """The cluster summary should pick the most-recently-deposited assembled isolate."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        r = conn.execute(
            "SELECT latest_assembly_json FROM cluster_summary WHERE pds_acc = ?",
            ("PDS_LIST_001.1",),
        ).fetchone()
    la = json.loads(r["latest_assembly_json"])
    assert la is not None
    # Among cluster 001's assembled isolates, PDT_HUMAN_001 has the most-recent
    # target_creation_date (2026-05-15). It has an asm_acc (GCA_001).
    assert la["asm_acc"] == "GCA_001"
    assert la["pdt"] == "PDT_HUMAN_001.1"
    assert la["target_creation_date"] == "2026-05-15"


def test_summarize_latest_assembly_skips_unassembled():
    """Isolates without an asm_acc must NOT be picked as the latest assembly."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    # In the fixture, PDT_OTHER_001.1 has the most-recent target_creation_date
    # in cluster 002 (2026-05-19). It has GCA_030 so SHOULD be the latest.
    with db.connect(db_path) as conn:
        r = conn.execute(
            "SELECT latest_assembly_json FROM cluster_summary WHERE pds_acc = ?",
            ("PDS_LIST_002.1",),
        ).fetchone()
    la = json.loads(r["latest_assembly_json"])
    assert la is not None
    assert la["asm_acc"] == "GCA_030"
    assert la["pdt"] == "PDT_OTHER_001.1"


def test_prose_ncbi_link_for():
    from src.render import prose
    # Prefer BioSample when available
    assert prose.ncbi_link_for("PDT_X", "SAMN123") == "https://www.ncbi.nlm.nih.gov/biosample/SAMN123"
    # Fall back to Pathogen Detection isolate browser otherwise
    assert prose.ncbi_link_for("PDT_X", None) == "https://www.ncbi.nlm.nih.gov/pathogens/isolates/#PDT_X"
    # Empty input → "#"
    assert prose.ncbi_link_for(None, None) == "#"


def test_summarize_histogram():
    """Per-year histogram should bin isolates correctly."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        r = conn.execute(
            "SELECT histogram_json, histogram_max_year_count "
            "FROM cluster_summary WHERE pds_acc = ?",
            ("PDS_LIST_001.1",),
        ).fetchone()
    hist = json.loads(r["histogram_json"])
    # PDS_LIST_001.1 dates: 10/15/25, 11/1/25 (x2), 12/10/25, 3/28/26, 4/10/26, 4/15/26
    # 2025 = 4 isolates (1 human + 3 nonhuman); 2026 = 3 isolates (all human)
    by_year = {h["year"]: h for h in hist}
    assert 2025 in by_year
    assert 2026 in by_year
    assert by_year[2025]["n_human"] == 1
    assert by_year[2025]["n_nonhuman"] == 3
    assert by_year[2026]["n_human"] == 3
    assert by_year[2026]["n_nonhuman"] == 0
    assert r["histogram_max_year_count"] == 4


def test_summarize_admin1():
    """Admin1 breakdown should preserve state-level specificity."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        r = conn.execute(
            "SELECT admin1_json FROM cluster_summary WHERE pds_acc = ?",
            ("PDS_LIST_001.1",),
        ).fetchone()
    admin1 = json.loads(r["admin1_json"])
    assert "USA" in admin1["by_country"]
    states = {a["admin1"]: a["n"] for a in admin1["by_country"]["USA"]}
    assert states.get("Texas") == 3
    assert states.get("Maryland") == 1
    assert states.get("California") == 1


def test_prose_format_geography():
    from src.render import prose
    assert prose.format_geography({}) == "—"
    data = {
        "by_country": {
            "USA": [
                {"admin1": "Texas", "n": 3},
                {"admin1": "Maryland", "n": 1},
            ],
        },
        "unspecified": {},
    }
    result = prose.format_geography(data)
    assert "USA" in result
    assert "Texas (3)" in result
    assert "Maryland (1)" in result
    data = {"by_country": {}, "unspecified": {"USA": 7}}
    result = prose.format_geography(data)
    assert "USA" in result
    assert "7" in result
    assert "location not specified" in result

    # Fix (1): redundant admin1 (where admin1 == country) should fold into unspecified
    data = {
        "by_country": {"United Kingdom": [{"admin1": "United Kingdom", "n": 5}]},
        "unspecified": {},
    }
    result = prose.format_geography(data)
    # "United Kingdom: United Kingdom (5)" should NOT appear; should be unspecified
    assert "United Kingdom: United Kingdom" not in result
    assert "location not specified" in result

    # Cap: more than max_countries should produce "+N more"
    data = {
        "by_country": {},
        "unspecified": {f"Country{i}": i for i in range(1, 12)},  # 11 countries
    }
    result = prose.format_geography(data)
    assert "+" in result
    assert "more countries" in result


def test_ifsac_excludes_human_isolates():
    """IFSAC counts should only reflect nonhuman isolates.

    A human isolate's IFSAC field (when present) describes reported
    food exposure, not isolate source — semantically different.
    """
    from src import signals
    members = [
        {"source_category": "Human", "ifsac_category": "Dairy"},        # should NOT count
        {"source_category": "Human", "ifsac_category": "Dairy"},        # should NOT count
        {"source_category": "Food",  "ifsac_category": "Meat-Poultry"}, # SHOULD count
        {"source_category": "Food",  "ifsac_category": "Meat-Poultry"}, # SHOULD count
    ]
    result = signals.compute_ifsac_summary(members)
    assert result["top"]["category"] == "Meat-Poultry"
    assert result["top"]["n"] == 2
    # Should NOT see Dairy (it's only on human isolates)
    cats = [c["category"] for c in result["categories"]]
    assert "Dairy" not in cats


def test_prose_format_oldest_isolate():
    from src.render import prose
    assert prose.format_oldest_isolate(None) == "—"
    assert prose.format_oldest_isolate({}) == "—"
    d = {
        "pdt": "PDT_HUMAN_001.1",
        "date_raw": "2026-04-15",
        "geo": "USA: Maryland",
        "source": "blood",
    }
    result = prose.format_oldest_isolate(d)
    assert "2026-04-15" in result
    assert "USA: Maryland" in result
    assert "blood" in result
    assert "PDT_HUMAN_001.1" in result


def test_histogram_svg_renders():
    """SVG histogram should render without error for typical data."""
    from src.render import histogram
    hist = [
        {"year": 2020, "n_human": 3, "n_nonhuman": 5},
        {"year": 2021, "n_human": 7, "n_nonhuman": 2},
        {"year": 2022, "n_human": 1, "n_nonhuman": 1},
    ]
    svg = histogram.render_histogram_svg(hist)
    assert svg.startswith("<svg")
    assert svg.endswith("</svg>")
    assert "rect" in svg
    assert "2020" in svg
    assert "2022" in svg


def test_histogram_svg_empty():
    """Empty histogram should render a graceful placeholder."""
    from src.render import histogram
    svg = histogram.render_histogram_svg([])
    assert svg.startswith("<svg")
    assert "No dated isolates" in svg


def test_histogram_svg_fills_year_gaps():
    """Years between min and max should appear even with zero counts."""
    from src.render import histogram
    hist = [
        {"year": 2020, "n_human": 1, "n_nonhuman": 0},
        {"year": 2023, "n_human": 2, "n_nonhuman": 0},
    ]
    svg = histogram.render_histogram_svg(hist)
    assert "2020" in svg
    assert "2023" in svg


# ---------------------------------------------------------------------------
# Helper inference tests
# ---------------------------------------------------------------------------

def test_infer_source_category_human():
    assert parse.infer_source_category("clinical", "Homo sapiens", None) == "Human"
    assert parse.infer_source_category(None, "Homo sapiens", None) == "Human"
    assert parse.infer_source_category(None, None, "blood") == "Human"
    assert parse.infer_source_category(None, None, "stool sample") == "Human"


def test_infer_source_category_food():
    assert parse.infer_source_category(None, None, "ground beef") == "Food"
    assert parse.infer_source_category(None, None, "deli meat - turkey") == "Food"
    assert parse.infer_source_category(None, None, "raw milk") == "Food"


def test_infer_source_category_animal():
    # Pure animal host with no food context → Animal
    assert parse.infer_source_category(None, "Cattle", None) == "Animal"
    assert parse.infer_source_category(None, "swine", "feces") == "Animal"  # animal host, non-food source
    # Animal host with food context → Food wins (it's a food signal)
    assert parse.infer_source_category(None, "Cattle", "milk") == "Food"
    assert parse.infer_source_category(None, "Chicken", "ground meat") == "Food"


def test_infer_source_category_environment():
    assert parse.infer_source_category(None, None, "wastewater") == "Environment"
    assert parse.infer_source_category(None, None, "soil sample") == "Environment"
    assert parse.infer_source_category(None, None, "drain swab") == "Environment"


def test_infer_source_category_unknown():
    assert parse.infer_source_category(None, None, None) == "Unknown"
    assert parse.infer_source_category(None, None, "some random thing") == "Unknown"


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
            print(f"  FAIL  {name}: {type(e).__name__}: {e}")
            failed += 1
    print()
    print(f"{passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
