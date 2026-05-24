"""Orchestrate the ingest pipeline.

For each configured pathogen:
  1. Fetch NCBI files (with cache freshness)
  2. Parse cluster_list → PDT→PDS map
  3. Parse metadata.tsv → IsolateRows
  4. Upsert into isolates table
  5. Parse amr.metadata.tsv → AmrRows
  6. Upsert into isolate_amr table
  7. Record the release in pdg_releases

Then once across all pathogens:
  8. Rebuild cluster_summary materialized view
"""

from __future__ import annotations

import argparse
import logging
import sys

from .. import config, db
from . import fetch, parse, summarize, upsert

log = logging.getLogger("pathogen-watch-ingest")


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def run() -> int:
    snapshots = fetch.fetch_all()
    if not snapshots:
        log.error("No snapshots fetched; aborting")
        return 1

    db.init_db()
    with db.connect() as conn:
        for snap in snapshots:
            log.info("=== %s (release %s) ===", snap.pathogen, snap.pdg_release)

            # Phase 1: cluster_list
            pds_map = parse.load_pds_map(snap.cluster_list_path)

            # Phase 2: metadata → isolates
            rows = parse.iter_isolates(
                snap.metadata_path,
                snap.pathogen,
                pds_map,
                snap.pdg_release,
            )
            n_iso = upsert.upsert_isolates(conn, rows)

            # Phase 3: AMR
            amr_rows = parse.iter_amr_rows(snap.amr_path)
            n_amr = upsert.upsert_amr(conn, amr_rows)

            # Phase 4: release record
            upsert.record_release(
                conn,
                pathogen=snap.pathogen,
                pdg_release=snap.pdg_release,
                metadata_url=snap.metadata_url,
                metadata_bytes=snap.metadata_bytes,
                cluster_list_url=snap.cluster_list_url,
                cluster_list_bytes=snap.cluster_list_bytes,
                amr_url=snap.amr_url,
                amr_bytes=snap.amr_bytes,
            )

            log.info("[%s] ingested: %d isolates, %d AMR gene calls",
                     snap.pathogen, n_iso, n_amr)

        # Phase 5: materialize cluster_summary
        summarize.materialize_cluster_summary(conn)

    log.info("Ingest complete.")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Pathogen Watch v2 ingest")
    ap.add_argument("-v", "--verbose", action="store_true", help="DEBUG logging")
    args = ap.parse_args(argv)
    setup_logging(args.verbose)
    return run()


if __name__ == "__main__":
    sys.exit(main())
