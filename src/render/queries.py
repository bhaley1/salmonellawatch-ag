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
               ct.mlst_representative_pdt, ct.mlst_error,
               ct.sistr_serovar, ct.sistr_antigenic_formula,
               ct.sistr_serogroup, ct.sistr_soc, ct.sistr_error,
               COALESCE(ct.consensus_serovar,
                   (SELECT serovar FROM isolates
                    WHERE pds_acc = cs.pds_acc
                    AND serovar IS NOT NULL AND serovar != ''
                    GROUP BY serovar ORDER BY COUNT(*) DESC LIMIT 1)
               ) AS effective_serovar
        FROM cluster_summary cs
        LEFT JOIN cluster_typing ct
               ON cs.pathogen = ct.pathogen AND cs.pds_acc = ct.pds_acc
        WHERE {where}
        ORDER BY cs.new_humans_in_window DESC, cs.latest_target_creation_date DESC
        LIMIT ?
    """, params).fetchall()

    # Pre-fetch AMR genotype summaries for active clusters
    pds_accs = [r["pds_acc"] for r in rows]
    placeholders = ",".join("?" * len(pds_accs))
    
    amr_by_cluster: dict[str, list[str]] = {}
    if pds_accs:
        amr_rows = conn.execute(f"""
            SELECT pds_acc, amr_genotypes, COUNT(*) as n
            FROM isolates
            WHERE pds_acc IN ({placeholders})
            AND amr_genotypes IS NOT NULL AND amr_genotypes != ''
            GROUP BY pds_acc, amr_genotypes
            ORDER BY pds_acc, n DESC
        """, pds_accs).fetchall()
        for ar in amr_rows:
            pds = ar["pds_acc"]
            if pds not in amr_by_cluster:
                amr_by_cluster[pds] = []
            amr_by_cluster[pds].append(ar["amr_genotypes"])

    # Pre-fetch recent non-human isolates (last 60 days) for active clusters
    nonhuman_by_cluster: dict[str, list[dict]] = {}
    if pds_accs:
        nh_rows = conn.execute(f"""
            SELECT pds_acc, epi_type, isolation_source, host,
                   geo_country, ifsac_category, target_creation_date,
                   COUNT(*) as n
            FROM isolates
            WHERE pds_acc IN ({placeholders})
            AND epi_type = 'environmental/other'
            AND target_creation_date >= date('now', '-60 days')
            GROUP BY pds_acc, geo_country, ifsac_category, isolation_source
            ORDER BY pds_acc, n DESC
        """, pds_accs).fetchall()
        for nr in nh_rows:
            pds = nr["pds_acc"]
            if pds not in nonhuman_by_cluster:
                nonhuman_by_cluster[pds] = []
            nonhuman_by_cluster[pds].append(dict(nr))

    # Pre-fetch source breakdown per cluster
    source_by_cluster: dict[str, dict] = {}
    if pds_accs:
        src_rows = conn.execute(f"""
            SELECT pds_acc, ifsac_category, COUNT(*) as n
            FROM isolates
            WHERE pds_acc IN ({placeholders})
            AND ifsac_category IS NOT NULL
            GROUP BY pds_acc, ifsac_category
        """, pds_accs).fetchall()
        for sr in src_rows:
            pds = sr["pds_acc"]
            if pds not in source_by_cluster:
                source_by_cluster[pds] = {}
            source_by_cluster[pds][sr["ifsac_category"]] = sr["n"]

    # Pre-fetch country-keyed human cases for filtering
    country_human_by_cluster: dict[str, dict] = {}
    if pds_accs:
        ch_rows = conn.execute(f"""
            SELECT pds_acc, geo_country,
                   pdt_acc, biosample_acc, collection_date,
                   collection_date_raw, geo_loc_name,
                   target_creation_date
            FROM isolates
            WHERE pds_acc IN ({placeholders})
            AND source_category = 'Human'
            AND collection_date >= date('now', '-60 days')
            ORDER BY pds_acc, collection_date DESC
        """, pds_accs).fetchall()
        for cr in ch_rows:
            pds = cr["pds_acc"]
            country = cr["geo_country"] or "Not Provided"
            if pds not in country_human_by_cluster:
                country_human_by_cluster[pds] = {}
            if country not in country_human_by_cluster[pds]:
                country_human_by_cluster[pds][country] = []
            country_human_by_cluster[pds][country].append({
                "pdt": cr["pdt_acc"],
                "biosample": cr["biosample_acc"],
                "collection_date": cr["collection_date_raw"] or cr["collection_date"],
                "geo": cr["geo_loc_name"],
                "geo_country": cr["geo_country"] or "",
            })

    # Pre-fetch country-keyed AMR genotypes
    country_amr_by_cluster: dict[str, dict] = {}
    if pds_accs:
        ca_rows = conn.execute(f"""
            SELECT pds_acc, geo_country, amr_genotypes, COUNT(*) as n
            FROM isolates
            WHERE pds_acc IN ({placeholders})
            AND amr_genotypes IS NOT NULL AND amr_genotypes != ''
            GROUP BY pds_acc, geo_country, amr_genotypes
            ORDER BY pds_acc, geo_country, n DESC
        """, pds_accs).fetchall()
        for ca in ca_rows:
            pds = ca["pds_acc"]
            country = ca["geo_country"] or "Not Provided"
            if pds not in country_amr_by_cluster:
                country_amr_by_cluster[pds] = {}
            if country not in country_amr_by_cluster[pds]:
                country_amr_by_cluster[pds][country] = []
            country_amr_by_cluster[pds][country].append(
                {"amr": ca["amr_genotypes"], "n": ca["n"]}
            )

    # Pre-fetch country-keyed source breakdown
    country_src_by_cluster: dict[str, dict] = {}
    if pds_accs:
        cs_rows = conn.execute(f"""
            SELECT pds_acc, geo_country, ifsac_category,
                   isolation_source, source_category, COUNT(*) as n
            FROM isolates
            WHERE pds_acc IN ({placeholders})
            GROUP BY pds_acc, geo_country, ifsac_category, source_category
            ORDER BY pds_acc, geo_country, n DESC
        """, pds_accs).fetchall()
        for cs in cs_rows:
            pds = cs["pds_acc"]
            country = cs["geo_country"] or "Not Provided"
            if pds not in country_src_by_cluster:
                country_src_by_cluster[pds] = {}
            if country not in country_src_by_cluster[pds]:
                country_src_by_cluster[pds][country] = []
            country_src_by_cluster[pds][country].append({
                "cat": cs["source_category"] or "Unknown",
                "src": cs["ifsac_category"] or cs["isolation_source"] or "(unspecified)",
                "n": cs["n"],
            })

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

        # SOC flag from SISTR
        # SOC: use SISTR result if available, else fall back to consensus_serovar
        EXACT_SOC = {
            "enteritidis","typhimurium","heidelberg","infantis","newport",
            "uganda","braenderup","muenchen","montevideo","javiana","reading",
            "dublin","oranienburg","potsdam","thompson","saintpaul","hadar",
            "schwarzengrund","anatum","berta"
        }
        SOC_FORMULAS = {"4:i:-","4,[5]:i:-","4,5:i:-","4,5,12:i:-","1,4,[5],12:i:-"}
        # Also match monophasic by serovar name patterns
        SOC_SEROVAR_PATTERNS = ["i 4", "i,4", "4,[5],12:i", "monophasic"]
        if d.get("sistr_soc") == 1:
            d["is_soc"] = True
        elif d.get("sistr_soc") == 0:
            # SISTR ran but returned not-SOC; double-check display serotype
            import re as _re3
            sv_raw = (d.get("effective_serovar") or d.get("consensus_serovar") or "")
            sv_check = _re3.sub(r"(?i)^salmonella\s+(enterica\s+subsp\.\s+enterica\s+serovar\s+)?", "", sv_raw).strip().lower()
            formula_check = (d.get("sistr_antigenic_formula") or "").lower()
            d["is_soc"] = (sv_check in EXACT_SOC or
                          any(f in formula_check for f in SOC_FORMULAS) or
                          any(p in sv_check for p in SOC_SEROVAR_PATTERNS))
        else:
            import re as _re2
            sv_raw = (d.get("effective_serovar") or d.get("consensus_serovar") or "")
            sv = _re2.sub(r"(?i)^salmonella\s+(enterica\s+subsp\.\s+enterica\s+serovar\s+)?", "", sv_raw).strip().lower()
            formula = (d.get("sistr_antigenic_formula") or "").lower()
            d["is_soc"] = (sv in EXACT_SOC or 
                          any(f in formula for f in SOC_FORMULAS) or
                          any(p in sv for p in SOC_SEROVAR_PATTERNS))

        # Attach pre-fetched aggregations
        pds = d["pds_acc"]
        d["amr_genotype_list"] = amr_by_cluster.get(pds, [])
        d["country_human_cases"] = country_human_by_cluster.get(pds, {})
        d["country_amr"] = country_amr_by_cluster.get(pds, {})
        d["country_src"] = country_src_by_cluster.get(pds, {})
        d["recent_nonhuman"] = nonhuman_by_cluster.get(pds, [])

        # Source breakdown bucketed for display
        raw_src = source_by_cluster.get(pds, {})
        src = {"beef": 0, "poultry": 0, "produce": 0, "rte": 0,
               "water": 0, "swine": 0, "other": {}}
        for cat, n in raw_src.items():
            cl = cat.lower()
            if any(x in cl for x in ["beef", "bovine", "cow"]):
                src["beef"] += n
            elif any(x in cl for x in ["chicken", "turkey", "poultry", "broiler"]):
                src["poultry"] += n
            elif any(x in cl for x in ["vegetable", "fruit", "produce", "leafy", "sprout", "herb", "nut"]):
                src["produce"] += n
            elif any(x in cl for x in ["ready-to-eat", "rte", "deli"]):
                src["rte"] += n
            elif any(x in cl for x in ["water", "aquatic", "stream"]):
                src["water"] += n
            elif any(x in cl for x in ["pork", "swine", "pig"]):
                src["swine"] += n
            elif "human" not in cl and "clinical" not in cl:
                src["other"][cat] = src["other"].get(cat, 0) + n
        d["source_breakdown"] = src
        d["sistr_formula"] = d.get("sistr_antigenic_formula")

        # Serotype: SISTR > effective_serovar (consensus or isolates fallback)
        import re as _re
        sistr_sv = d.get("sistr_serovar")
        if sistr_sv and sistr_sv not in ("", "-"):
            d["consensus_serotype"] = sistr_sv
        else:
            raw = d.get("effective_serovar") or d.get("consensus_serovar") or ""
            clean = _re.sub(r"(?i)^salmonella\s+(enterica\s+subsp\.\s+enterica\s+serovar\s+)?", "", raw).strip()
            if clean and not _re.match(r"^[A-Z]", clean):
                clean = clean.capitalize()
            d["consensus_serotype"] = clean if clean else None
        d["consensus_serotype_n"] = d.get("consensus_serovar_n")
        d["consensus_serotype_total"] = d.get("consensus_serovar_total")
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
