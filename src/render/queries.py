"""SQL queries used by the renderer.

Centralizing queries here keeps the templates clean and makes the data
contract between ingest and render explicit. All queries return plain
dicts ready for Jinja consumption.
"""

from __future__ import annotations

import json
import sqlite3
from typing import Any


def get_recent_activity_clusters(
    conn: sqlite3.Connection,
    pathogen: str | None = None,
    limit: int = 500,
) -> list[dict[str, Any]]:
    """Clusters with new human additions in the last 30 days.

    Returns the primary surveillance data for the dashboard, ordered by
    number of new humans descending. Each row includes the per-cluster
    summary plus the typed/decoded JSON payloads, and (where available)
    the cluster's consensus serovar and MLST sequence type.
    """
    where = "cs.new_humans_in_window > 0"
    params: list[Any] = []
    if pathogen:
        where += " AND cs.pathogen = ?"
        params.append(pathogen)
    params.append(limit)

    rows = conn.execute(f"""
        SELECT cs.pathogen, cs.pds_acc,
               cs.n_total, cs.n_human, cs.n_nonhuman,
               cs.n_food, cs.n_animal, cs.n_environment,
               cs.earliest_collection_date, cs.latest_collection_date,
               cs.earliest_target_creation_date, cs.latest_target_creation_date,
               cs.temporal_span_days,
               cs.oldest_isolate_json, cs.oldest_human_json, cs.oldest_nonhuman_json,
               cs.histogram_json, cs.histogram_max_year_count,
               cs.deposit_lag_median, cs.deposit_lag_mean,
               cs.host_summary_json,
               cs.map_locations_json, cs.latest_assembly_json,
               cs.signals_json,
               cs.countries_json, cs.admin1_json, cs.source_summary_json,
               cs.new_humans_in_window,
               cs.new_humans_30d,
               cs.new_humans_15d,
               cs.new_humans_in_window_pdts_json,
               cs.new_humans_in_window_dates_json,
               cs.window_days,
               ct.consensus_serovar,
               ct.consensus_serovar_n,
               ct.consensus_serovar_total,
               ct.mlst_scheme, ct.mlst_st, ct.mlst_alleles,
               ct.mlst_representative_pdt, ct.mlst_error
        FROM cluster_summary cs
        LEFT JOIN cluster_typing ct
               ON cs.pathogen = ct.pathogen AND cs.pds_acc = ct.pds_acc
        WHERE {where}
        ORDER BY cs.new_humans_in_window DESC, cs.latest_target_creation_date DESC
        LIMIT ?
    """, params).fetchall()

    out: list[dict[str, Any]] = []
    for r in rows:
        d = dict(r)
        d["countries"] = json.loads(d.pop("countries_json") or "[]")
        d["admin1"] = json.loads(d.pop("admin1_json") or '{"by_country":{},"unspecified":{}}')
        d["source_summary"] = json.loads(d.pop("source_summary_json") or "[]")
        d["new_humans_pdts"] = json.loads(d.pop("new_humans_in_window_pdts_json") or "[]")
        d["new_humans_dates"] = json.loads(d.pop("new_humans_in_window_dates_json") or "[]")
        d["oldest_any"] = json.loads(d.pop("oldest_isolate_json") or "null")
        d["oldest_human"] = json.loads(d.pop("oldest_human_json") or "null")
        d["oldest_nonhuman"] = json.loads(d.pop("oldest_nonhuman_json") or "null")
        d["histogram"] = json.loads(d.pop("histogram_json") or "[]")
        d["map_locations"] = json.loads(d.pop("map_locations_json") or "[]")
        d["latest_assembly"] = json.loads(d.pop("latest_assembly_json") or "null")
        d["host_summary"] = json.loads(d.pop("host_summary_json") or "[]")
        d["signals"] = json.loads(d.pop("signals_json") or "{}")
        d["mlst_alleles_dict"] = json.loads(d["mlst_alleles"]) if d.get("mlst_alleles") else {}
        out.append(d)
    return out


def get_totals(conn: sqlite3.Connection) -> dict[str, Any]:
    """Site-wide totals shown in tiles."""
    row = conn.execute("""
        SELECT
            COUNT(*) AS n_clusters,
            SUM(CASE WHEN n_human > 0 THEN 1 ELSE 0 END) AS n_human_clusters,
            SUM(CASE WHEN n_human > 0 AND n_nonhuman > 0 THEN 1 ELSE 0 END) AS n_mixed,
            SUM(CASE WHEN new_humans_in_window > 0 THEN 1 ELSE 0 END) AS n_active,
            SUM(new_humans_in_window) AS n_new_humans_window
        FROM cluster_summary
    """).fetchone()
    return dict(row) if row else {}


def get_pathogen_counts(conn: sqlite3.Connection) -> list[dict[str, Any]]:
    """For chip filters: count of active clusters per pathogen."""
    rows = conn.execute("""
        SELECT pathogen,
               COUNT(*) AS n_clusters,
               SUM(CASE WHEN new_humans_in_window > 0 THEN 1 ELSE 0 END) AS n_active
        FROM cluster_summary
        GROUP BY pathogen
        ORDER BY pathogen
    """).fetchall()
    return [dict(r) for r in rows]


def get_latest_release(conn: sqlite3.Connection) -> dict[str, Any] | None:
    """Most recent PDG release we've ingested across all pathogens."""
    row = conn.execute("""
        SELECT pathogen, pdg_release, ingested_at
        FROM pdg_releases
        ORDER BY ingested_at DESC
        LIMIT 1
    """).fetchone()
    return dict(row) if row else None
