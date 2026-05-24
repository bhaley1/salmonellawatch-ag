"""Orchestrate the cluster-level subtyping pipeline.

Two operations:

  1. Consensus serovar — fast, pure DB aggregation. Runs against all
     clusters every time; result is stored in cluster_typing.

  2. MLST typing — slow (assembly download + BLAST). Runs only against
     clusters that don't yet have a typing record. Budget-capped so
     first runs don't exhaust GitHub Actions' 6-hour limit; convergence
     happens over ~1 week of daily runs.

Both operations write to the same cluster_typing table via upsert.
"""

from __future__ import annotations

import argparse
import logging
import sqlite3
import sys
from datetime import datetime

from .. import config, db
from . import assembly, mlst, representative, serovar

log = logging.getLogger("pathogen-watch-subtype")


# Cap typings per run so initial backfills don't exceed Actions timeouts.
# At ~30s per cluster (mostly assembly download + BLAST), 500 = ~4 hours.
DEFAULT_MLST_BUDGET = 500


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def _upsert_typing(
    conn: sqlite3.Connection,
    pathogen: str,
    pds_acc: str,
    serovar_consensus: serovar.SerovarConsensus,
    mlst_result: mlst.MLSTResult | None,
    representative_pdt: str | None,
) -> None:
    """Upsert a cluster_typing row.

    Serovar and MLST come from different passes. We want each pass to
    update only the fields it has data for, never wiping out the other.
    The serovar pass passes mlst_result=None; the MLST pass passes a
    placeholder SerovarConsensus(None, 0, 0) which we detect via
    n_total_with_serovar==0.
    """
    now = datetime.utcnow().isoformat(timespec="seconds")
    has_serovar_data = serovar_consensus.n_total_with_serovar > 0
    has_mlst_data = mlst_result is not None

    # Read existing row, if any, so we can preserve the other side's data
    existing = conn.execute(
        "SELECT * FROM cluster_typing WHERE pathogen = ? AND pds_acc = ?",
        (pathogen, pds_acc),
    ).fetchone()

    # Resolve final values: prefer new data when this pass has it,
    # otherwise fall back to whatever the existing row has.
    if has_serovar_data:
        final_serovar = serovar_consensus.consensus_serovar
        final_serovar_n = serovar_consensus.n_agreed
        final_serovar_total = serovar_consensus.n_total_with_serovar
    else:
        final_serovar = existing["consensus_serovar"] if existing else None
        final_serovar_n = existing["consensus_serovar_n"] if existing else 0
        final_serovar_total = existing["consensus_serovar_total"] if existing else 0

    if has_mlst_data:
        final_mlst_scheme = mlst_result.scheme
        final_mlst_st = mlst_result.st
        final_mlst_alleles = mlst.alleles_to_json(mlst_result.alleles) if mlst_result.alleles else None
        final_mlst_rep = representative_pdt
        final_mlst_error = mlst_result.error
    else:
        final_mlst_scheme = existing["mlst_scheme"] if existing else None
        final_mlst_st = existing["mlst_st"] if existing else None
        final_mlst_alleles = existing["mlst_alleles"] if existing else None
        final_mlst_rep = existing["mlst_representative_pdt"] if existing else None
        final_mlst_error = existing["mlst_error"] if existing else None

    conn.execute("""
        INSERT INTO cluster_typing (
            pathogen, pds_acc,
            consensus_serovar, consensus_serovar_n, consensus_serovar_total,
            mlst_scheme, mlst_st, mlst_alleles,
            mlst_representative_pdt, mlst_error,
            typed_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(pathogen, pds_acc) DO UPDATE SET
            consensus_serovar       = excluded.consensus_serovar,
            consensus_serovar_n     = excluded.consensus_serovar_n,
            consensus_serovar_total = excluded.consensus_serovar_total,
            mlst_scheme             = excluded.mlst_scheme,
            mlst_st                 = excluded.mlst_st,
            mlst_alleles            = excluded.mlst_alleles,
            mlst_representative_pdt = excluded.mlst_representative_pdt,
            mlst_error              = excluded.mlst_error,
            typed_at                = excluded.typed_at
    """, (
        pathogen, pds_acc,
        final_serovar, final_serovar_n, final_serovar_total,
        final_mlst_scheme, final_mlst_st, final_mlst_alleles,
        final_mlst_rep, final_mlst_error,
        now,
    ))


def run_serovar_consensus(conn: sqlite3.Connection) -> int:
    """Update consensus_serovar for every cluster. Fast (pure SQL aggregation)."""
    log.info("Computing consensus serovars for all clusters...")
    results = serovar.compute_all_cluster_serovars(conn)

    n_with_consensus = 0
    for pathogen, pds_acc, consensus in results:
        if consensus.consensus_serovar:
            n_with_consensus += 1
        # Upsert preserves any existing MLST data for this cluster
        _upsert_typing(
            conn,
            pathogen=pathogen,
            pds_acc=pds_acc,
            serovar_consensus=consensus,
            mlst_result=None,
            representative_pdt=None,
        )
    conn.commit()
    log.info(
        "Serovar consensus: %d clusters typed (%d with a non-null consensus)",
        len(results), n_with_consensus,
    )
    return len(results)


def run_mlst_backfill(
    conn: sqlite3.Connection,
    budget: int = DEFAULT_MLST_BUDGET,
    pathogen_filter: str | None = None,
) -> dict[str, int]:
    """Type up to `budget` clusters that don't yet have MLST results.

    Returns counters by outcome: {'typed', 'no_assembly', 'tool_missing',
    'fetch_failed', 'mlst_error'}.
    """
    counters = {
        "typed": 0,
        "no_assembly": 0,
        "tool_missing": 0,
        "fetch_failed": 0,
        "mlst_error": 0,
    }

    if not mlst.is_mlst_installed():
        log.warning(
            "mlst CLI not installed. MLST backfill skipping; "
            "install Torsten Seemann's mlst tool to enable typing."
        )
        counters["tool_missing"] = -1  # sentinel: not even attempted
        return counters

    pending = representative.list_clusters_needing_typing(
        conn, pathogen=pathogen_filter, limit=budget,
    )
    log.info("MLST backfill: %d clusters in queue (budget=%d)", len(pending), budget)

    for pathogen, pds_acc in pending:
        # 1. Pick a representative
        rep = representative.pick_representative(conn, pathogen, pds_acc)
        if not rep:
            log.info("[%s %s] no assembly available; marking no_assembly", pathogen, pds_acc)
            _upsert_typing(
                conn, pathogen, pds_acc,
                serovar_consensus=serovar.SerovarConsensus(None, 0, 0),
                mlst_result=mlst.MLSTResult(
                    scheme=mlst.PATHOGEN_SCHEME.get(pathogen),
                    st=None, alleles={}, error="no_assembly",
                ),
                representative_pdt=None,
            )
            counters["no_assembly"] += 1
            conn.commit()
            continue

        scheme = mlst.PATHOGEN_SCHEME.get(pathogen)
        if not scheme:
            log.warning("No MLST scheme configured for %s", pathogen)
            continue

        # 2. Fetch the assembly
        fasta = assembly.fetch_assembly_fasta(rep.asm_acc)
        if not fasta:
            log.info("[%s %s] assembly fetch failed for %s", pathogen, pds_acc, rep.asm_acc)
            _upsert_typing(
                conn, pathogen, pds_acc,
                serovar_consensus=serovar.SerovarConsensus(None, 0, 0),
                mlst_result=mlst.MLSTResult(
                    scheme=scheme, st=None, alleles={},
                    error=f"fetch_failed:{rep.asm_acc}",
                ),
                representative_pdt=rep.pdt_acc,
            )
            counters["fetch_failed"] += 1
            conn.commit()
            continue

        # 3. Run MLST
        result = mlst.run_mlst(fasta, scheme)
        _upsert_typing(
            conn, pathogen, pds_acc,
            serovar_consensus=serovar.SerovarConsensus(None, 0, 0),
            mlst_result=result,
            representative_pdt=rep.pdt_acc,
        )
        conn.commit()

        if result.error == "tool_not_installed":
            counters["tool_missing"] += 1
            # Don't keep trying every cluster if the tool isn't there
            log.warning("Aborting MLST backfill: tool not installed")
            break
        elif result.st and not result.error:
            log.info("[%s %s] typed: %s (rep: %s)",
                     pathogen, pds_acc, result.st, rep.asm_acc)
            counters["typed"] += 1
        else:
            log.info("[%s %s] mlst error: %s", pathogen, pds_acc, result.error)
            counters["mlst_error"] += 1

    log.info("MLST backfill complete: %s", counters)
    return counters


def run(
    *,
    mlst_budget: int = DEFAULT_MLST_BUDGET,
    skip_mlst: bool = False,
    pathogen: str | None = None,
) -> int:
    db.init_db()
    with db.connect() as conn:
        run_serovar_consensus(conn)
        if not skip_mlst:
            run_mlst_backfill(conn, budget=mlst_budget, pathogen_filter=pathogen)
        else:
            log.info("MLST backfill skipped (--skip-mlst)")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Pathogen Watch v2 cluster subtyping")
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("--skip-mlst", action="store_true",
                    help="Only run consensus serovar; skip MLST")
    ap.add_argument("--mlst-budget", type=int, default=DEFAULT_MLST_BUDGET,
                    help="Cap number of MLST runs this invocation (default: 500)")
    ap.add_argument("--pathogen", help="Limit to one pathogen")
    args = ap.parse_args(argv)
    setup_logging(args.verbose)
    return run(
        mlst_budget=args.mlst_budget,
        skip_mlst=args.skip_mlst,
        pathogen=args.pathogen,
    )


if __name__ == "__main__":
    sys.exit(main())
