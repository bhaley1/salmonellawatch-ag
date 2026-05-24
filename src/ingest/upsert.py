"""Upsert parsed records into SQLite.

Uses INSERT...ON CONFLICT DO UPDATE so re-running the ingest pipeline is
idempotent and incremental: existing rows get their values refreshed,
new rows get inserted. Batches of 1000 for performance.
"""

from __future__ import annotations

import logging
import sqlite3
from datetime import datetime
from typing import Iterable

from .parse import AmrRow, IsolateRow

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Isolates
# ---------------------------------------------------------------------------
_ISOLATE_UPSERT = """
INSERT INTO isolates (
    pdt_acc, pathogen, pds_acc, epi_type, source_category,
    host, isolation_source, geo_loc_name, geo_country, geo_admin1,
    collection_date, collection_date_raw, target_creation_date,
    scientific_name, serovar, biosample_acc, asm_acc, sra_acc,
    asm_level, asm_stats_contig_n50,
    food_origin, ifsac_category, host_disease, bioproject_acc,
    pdg_release, last_seen_at
) VALUES (
    :pdt_acc, :pathogen, :pds_acc, :epi_type, :source_category,
    :host, :isolation_source, :geo_loc_name, :geo_country, :geo_admin1,
    :collection_date, :collection_date_raw, :target_creation_date,
    :scientific_name, :serovar, :biosample_acc, :asm_acc, :sra_acc,
    :asm_level, :asm_stats_contig_n50,
    :food_origin, :ifsac_category, :host_disease, :bioproject_acc,
    :pdg_release, :last_seen_at
)
ON CONFLICT(pdt_acc) DO UPDATE SET
    pds_acc              = excluded.pds_acc,
    epi_type             = excluded.epi_type,
    source_category      = excluded.source_category,
    host                 = excluded.host,
    isolation_source     = excluded.isolation_source,
    geo_loc_name         = excluded.geo_loc_name,
    geo_country          = excluded.geo_country,
    geo_admin1           = excluded.geo_admin1,
    collection_date      = excluded.collection_date,
    collection_date_raw  = excluded.collection_date_raw,
    target_creation_date = excluded.target_creation_date,
    scientific_name      = excluded.scientific_name,
    serovar              = excluded.serovar,
    biosample_acc        = excluded.biosample_acc,
    asm_acc              = excluded.asm_acc,
    sra_acc              = excluded.sra_acc,
    asm_level            = excluded.asm_level,
    asm_stats_contig_n50 = excluded.asm_stats_contig_n50,
    food_origin          = excluded.food_origin,
    ifsac_category       = excluded.ifsac_category,
    host_disease         = excluded.host_disease,
    bioproject_acc       = excluded.bioproject_acc,
    pdg_release          = excluded.pdg_release,
    last_seen_at         = excluded.last_seen_at
"""


def upsert_isolates(
    conn: sqlite3.Connection,
    rows: Iterable[IsolateRow],
    batch_size: int = 1000,
) -> int:
    """Upsert isolate rows in batches. Returns count inserted/updated."""
    now = datetime.utcnow().isoformat(timespec="seconds")
    batch: list[dict] = []
    count = 0

    def _flush():
        nonlocal batch
        if batch:
            conn.executemany(_ISOLATE_UPSERT, batch)
            conn.commit()
            batch = []

    for row in rows:
        batch.append({
            "pdt_acc":              row.pdt_acc,
            "pathogen":             row.pathogen,
            "pds_acc":              row.pds_acc,
            "epi_type":             row.epi_type,
            "source_category":      row.source_category,
            "host":                 row.host,
            "isolation_source":     row.isolation_source,
            "geo_loc_name":         row.geo_loc_name,
            "geo_country":          row.geo_country,
            "geo_admin1":           row.geo_admin1,
            "collection_date":      row.collection_date.isoformat() if row.collection_date else None,
            "collection_date_raw":  row.collection_date_raw,
            "target_creation_date": row.target_creation_date.isoformat() if row.target_creation_date else None,
            "scientific_name":      row.scientific_name,
            "serovar":              row.serovar,
            "biosample_acc":        row.biosample_acc,
            "asm_acc":              row.asm_acc,
            "sra_acc":              row.sra_acc,
            "asm_level":            row.asm_level,
            "asm_stats_contig_n50": row.asm_stats_contig_n50,
            "food_origin":          row.food_origin,
            "ifsac_category":       row.ifsac_category,
            "host_disease":         row.host_disease,
            "bioproject_acc":       row.bioproject_acc,
            "pdg_release":          row.pdg_release,
            "last_seen_at":         now,
        })
        count += 1
        if len(batch) >= batch_size:
            _flush()
    _flush()
    log.info("Upserted %d isolate rows", count)
    return count


# ---------------------------------------------------------------------------
# AMR
# ---------------------------------------------------------------------------
_AMR_UPSERT = """
INSERT INTO isolate_amr (pdt_acc, gene_symbol, element_subtype)
VALUES (:pdt_acc, :gene_symbol, :element_subtype)
ON CONFLICT(pdt_acc, gene_symbol) DO UPDATE SET
    element_subtype = excluded.element_subtype
"""


def upsert_amr(
    conn: sqlite3.Connection,
    rows: Iterable[AmrRow],
    batch_size: int = 5000,
) -> int:
    """Upsert AMR gene-call rows. AMR is row-per-gene so batches can be larger."""
    batch: list[dict] = []
    count = 0

    def _flush():
        nonlocal batch
        if batch:
            conn.executemany(_AMR_UPSERT, batch)
            conn.commit()
            batch = []

    for row in rows:
        batch.append({
            "pdt_acc":         row.pdt_acc,
            "gene_symbol":     row.gene_symbol,
            "element_subtype": row.element_subtype,
        })
        count += 1
        if len(batch) >= batch_size:
            _flush()
    _flush()
    log.info("Upserted %d AMR rows", count)
    return count


# ---------------------------------------------------------------------------
# Release log
# ---------------------------------------------------------------------------
def record_release(
    conn: sqlite3.Connection,
    pathogen: str,
    pdg_release: str,
    metadata_url: str,
    metadata_bytes: int,
    cluster_list_url: str | None,
    cluster_list_bytes: int | None,
    amr_url: str | None,
    amr_bytes: int | None,
) -> None:
    now = datetime.utcnow().isoformat(timespec="seconds")
    conn.execute(
        """
        INSERT INTO pdg_releases (
            pathogen, pdg_release, metadata_url, metadata_bytes,
            cluster_list_url, cluster_list_bytes, amr_url, amr_bytes,
            ingested_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(pathogen, pdg_release) DO UPDATE SET
            metadata_url       = excluded.metadata_url,
            metadata_bytes     = excluded.metadata_bytes,
            cluster_list_url   = excluded.cluster_list_url,
            cluster_list_bytes = excluded.cluster_list_bytes,
            amr_url            = excluded.amr_url,
            amr_bytes          = excluded.amr_bytes,
            ingested_at        = excluded.ingested_at
        """,
        (pathogen, pdg_release, metadata_url, metadata_bytes,
         cluster_list_url, cluster_list_bytes, amr_url, amr_bytes, now),
    )
    conn.commit()
