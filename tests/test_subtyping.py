"""Phase 2 tests: cluster-level subtyping.

These tests cover:
  - Consensus serovar computation (pure logic)
  - MLST output parsing (pure logic, no CLI required)
  - Representative selection (DB query logic)
  - The orchestrator end-to-end with mocked MLST and assembly fetching

The real `mlst` CLI is NOT required to run these tests.
"""

import json
import sqlite3
import sys
import tempfile
from datetime import date
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import config, db
from src.subtyping import mlst, representative, run as subtyping_run, serovar

# Reuse the Phase 1 fixture loader
sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_pipeline import fresh_db, load_fixture_into_db


# ---------------------------------------------------------------------------
# Consensus serovar
# ---------------------------------------------------------------------------

def test_serovar_consensus_all_agree():
    result = serovar.compute_consensus_serovar(["1/2a", "1/2a", "1/2a"])
    assert result.consensus_serovar == "1/2a"
    assert result.n_agreed == 3
    assert result.n_total_with_serovar == 3
    assert result.agreement_fraction == 1.0


def test_serovar_consensus_majority():
    result = serovar.compute_consensus_serovar(["1/2a", "1/2a", "4b"])
    assert result.consensus_serovar == "1/2a"
    assert result.n_agreed == 2
    assert result.n_total_with_serovar == 3
    assert abs(result.agreement_fraction - 2/3) < 0.01


def test_serovar_consensus_with_nulls():
    """Null/empty values are ignored, not counted as agreement."""
    result = serovar.compute_consensus_serovar([None, "", "1/2a", "1/2a"])
    assert result.consensus_serovar == "1/2a"
    assert result.n_agreed == 2
    assert result.n_total_with_serovar == 2


def test_serovar_consensus_all_null():
    result = serovar.compute_consensus_serovar([None, None, ""])
    assert result.consensus_serovar is None
    assert result.n_agreed == 0
    assert result.n_total_with_serovar == 0


def test_serovar_consensus_from_db():
    """End-to-end against the synthetic Listeria fixture."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        results = serovar.compute_all_cluster_serovars(conn)
    by_cluster = {(p, pds): c for p, pds, c in results}

    # PDS_LIST_001.1 has all members with serovar 1/2a
    c1 = by_cluster[("Listeria", "PDS_LIST_001.1")]
    assert c1.consensus_serovar == "1/2a"
    assert c1.n_agreed == 7
    assert c1.n_total_with_serovar == 7

    # PDS_LIST_002.1 has 2 members, both serovar 4b
    c2 = by_cluster[("Listeria", "PDS_LIST_002.1")]
    assert c2.consensus_serovar == "4b"
    assert c2.n_agreed == 2


# ---------------------------------------------------------------------------
# MLST output parsing
# ---------------------------------------------------------------------------

def test_mlst_parse_clean_st():
    output = "assembly.fna\tlmonocytogenes_2\t6\tabcZ(3)\tbglA(1)\tcat(1)\tdapE(1)\tdat(3)\tldh(1)\tlhkA(1)\n"
    result = mlst.parse_mlst_output(output, "lmonocytogenes_2")
    assert result.scheme == "lmonocytogenes_2"
    assert result.st == "ST6"
    assert result.alleles == {
        "abcZ": "3", "bglA": "1", "cat": "1",
        "dapE": "1", "dat": "3", "ldh": "1", "lhkA": "1",
    }
    assert result.error is None


def test_mlst_parse_untypeable():
    output = "assembly.fna\t-\t-\n"
    result = mlst.parse_mlst_output(output, "lmonocytogenes_2")
    assert result.st == "untypeable"
    assert result.error == "no_match"


def test_mlst_parse_novel():
    """Partial allele (~) means a novel variant; ST is flagged."""
    output = "assembly.fna\tlmonocytogenes_2\t6\tabcZ(~3)\tbglA(1)\tcat(1)\tdapE(1)\tdat(3)\tldh(1)\tlhkA(1)\n"
    result = mlst.parse_mlst_output(output, "lmonocytogenes_2")
    assert "novel" in result.st.lower()
    assert result.alleles["abcZ"] == "~3"


def test_mlst_parse_empty():
    result = mlst.parse_mlst_output("", "lmonocytogenes_2")
    assert result.st is None
    assert result.error == "empty_output"


def test_mlst_parse_malformed():
    result = mlst.parse_mlst_output("not a tsv\n", "lmonocytogenes_2")
    assert result.st is None
    assert result.error and result.error.startswith("malformed_output")


def test_mlst_alleles_to_json_stable():
    """JSON serialization should be deterministic for caching."""
    a = {"bglA": "1", "abcZ": "3"}
    b = {"abcZ": "3", "bglA": "1"}
    assert mlst.alleles_to_json(a) == mlst.alleles_to_json(b)


# ---------------------------------------------------------------------------
# Representative selection
# ---------------------------------------------------------------------------

def test_representative_picks_best_assembly_level():
    """Complete Genome > Chromosome > Contig, with n50 as tiebreaker."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        rep = representative.pick_representative(conn, "Listeria", "PDS_LIST_001.1")
    assert rep is not None
    # Three isolates have asm_acc + 'Complete Genome': HUMAN_001 (n50=2.9M),
    # FOOD_001 (n50=2.95M), and FOOD_002 (Contig, n50=450k). The highest n50
    # Complete Genome wins.
    assert rep.asm_level == "Complete Genome"
    assert rep.pdt_acc == "PDT_FOOD_001.1"  # n50=2,950,000, highest among Complete Genomes
    assert rep.asm_acc == "GCA_010"


def test_representative_skips_no_assembly():
    """Isolates without asm_acc should be skipped during selection."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        rep = representative.pick_representative(conn, "Listeria", "PDS_LIST_001.1")
    # The picked representative must have a non-null assembly
    assert rep.asm_acc and rep.asm_acc.strip()


def test_representative_returns_none_when_no_assemblies():
    """If no member of a cluster has an assembly, return None."""
    db_path = fresh_db()
    with db.connect(db_path) as conn:
        # Insert one isolate in a fake cluster with no assembly
        conn.execute(
            "INSERT INTO isolates (pdt_acc, pathogen, pds_acc, asm_acc, last_seen_at, pdg_release) "
            "VALUES ('PDT_X.1', 'Listeria', 'PDS_NOASM.1', NULL, '2026-05-21', 'X.1')",
        )
        conn.commit()
        rep = representative.pick_representative(conn, "Listeria", "PDS_NOASM.1")
    assert rep is None


def test_list_clusters_needing_typing():
    """Clusters with no cluster_typing row should appear in the queue."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        pending = representative.list_clusters_needing_typing(conn)
    # Both fixture clusters should be pending typing
    assert len(pending) == 2
    pds_accs = {p[1] for p in pending}
    assert "PDS_LIST_001.1" in pds_accs
    assert "PDS_LIST_002.1" in pds_accs


def test_list_clusters_needing_typing_skips_typed():
    """If a cluster already has mlst_st, it shouldn't reappear."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        conn.execute("""
            INSERT INTO cluster_typing (
                pathogen, pds_acc, mlst_scheme, mlst_st, typed_at
            ) VALUES ('Listeria', 'PDS_LIST_001.1', 'lmonocytogenes_2', 'ST6', '2026-05-21')
        """)
        conn.commit()
        pending = representative.list_clusters_needing_typing(conn)
    pds_accs = {p[1] for p in pending}
    assert "PDS_LIST_001.1" not in pds_accs
    assert "PDS_LIST_002.1" in pds_accs


def test_list_clusters_needing_typing_skips_errored():
    """Clusters with mlst_error set are NOT retried (avoid infinite retry)."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)
    with db.connect(db_path) as conn:
        conn.execute("""
            INSERT INTO cluster_typing (
                pathogen, pds_acc, mlst_error, typed_at
            ) VALUES ('Listeria', 'PDS_LIST_001.1', 'no_assembly', '2026-05-21')
        """)
        conn.commit()
        pending = representative.list_clusters_needing_typing(conn)
    pds_accs = {p[1] for p in pending}
    assert "PDS_LIST_001.1" not in pds_accs


# ---------------------------------------------------------------------------
# Orchestrator with mocked MLST
# ---------------------------------------------------------------------------

def test_orchestrator_with_mocked_mlst():
    """Run the full subtyping pipeline with mocked external calls."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)

    fake_fasta = Path("/tmp/fake.fna")

    def fake_fetch(asm_acc, cache_dir=None):
        return fake_fasta  # pretend any assembly downloads cleanly

    def fake_is_installed():
        return True

    def fake_run_mlst(fasta_path, scheme, timeout=300):
        return mlst.MLSTResult(
            scheme=scheme, st="ST6",
            alleles={"abcZ": "3", "bglA": "1"},
            error=None,
        )

    # Patch config.DB_PATH to point at our test DB
    with patch.object(config, "DB_PATH", db_path), \
         patch("src.subtyping.run.assembly.fetch_assembly_fasta", side_effect=fake_fetch), \
         patch("src.subtyping.run.mlst.is_mlst_installed", side_effect=fake_is_installed), \
         patch("src.subtyping.run.mlst.run_mlst", side_effect=fake_run_mlst):
        rc = subtyping_run.run(skip_mlst=False, mlst_budget=10)
    assert rc == 0

    with db.connect(db_path) as conn:
        rows = conn.execute("""
            SELECT pds_acc, consensus_serovar, mlst_st, mlst_error
            FROM cluster_typing ORDER BY pds_acc
        """).fetchall()
    assert len(rows) == 2
    by_pds = {r["pds_acc"]: dict(r) for r in rows}
    assert by_pds["PDS_LIST_001.1"]["consensus_serovar"] == "1/2a"
    assert by_pds["PDS_LIST_001.1"]["mlst_st"] == "ST6"
    assert by_pds["PDS_LIST_001.1"]["mlst_error"] is None
    assert by_pds["PDS_LIST_002.1"]["consensus_serovar"] == "4b"


def test_orchestrator_handles_missing_mlst_tool():
    """When mlst is not installed, serovar still runs; MLST records the error."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)

    def fake_not_installed():
        return False

    with patch.object(config, "DB_PATH", db_path), \
         patch("src.subtyping.run.mlst.is_mlst_installed", side_effect=fake_not_installed):
        rc = subtyping_run.run(skip_mlst=False, mlst_budget=10)
    assert rc == 0

    with db.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT pds_acc, consensus_serovar FROM cluster_typing"
        ).fetchall()
    # Serovars get populated even without the MLST tool
    assert len(rows) == 2
    assert all(r["consensus_serovar"] for r in rows)


def test_orchestrator_skip_mlst_flag():
    """--skip-mlst leaves cluster_typing populated with serovar only."""
    db_path = fresh_db()
    load_fixture_into_db(db_path)

    with patch.object(config, "DB_PATH", db_path):
        rc = subtyping_run.run(skip_mlst=True)
    assert rc == 0

    with db.connect(db_path) as conn:
        rows = conn.execute(
            "SELECT pds_acc, consensus_serovar, mlst_st FROM cluster_typing"
        ).fetchall()
    assert len(rows) == 2
    assert all(r["consensus_serovar"] for r in rows)
    assert all(r["mlst_st"] is None for r in rows)


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
