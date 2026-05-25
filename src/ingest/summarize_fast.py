"""Fast cluster_summary materialization using bulk SQL.

Replaces the slow per-cluster Python loop with bulk SQL aggregations.
Only the fields that truly require Python logic (signals, histogram,
map locations) are computed in a second pass over pre-fetched data.
"""

from __future__ import annotations

import json
import logging
import sqlite3
from datetime import date, datetime, timedelta

from .. import config

log = logging.getLogger(__name__)


def materialize_cluster_summary(
    conn: sqlite3.Connection,
    window_days: int | None = None,
    today: date | None = None,
) -> int:
    window_days = window_days or config.RECENT_WINDOW_DAYS
    today = today or date.today()
    cutoff_60 = (today - timedelta(days=window_days)).isoformat()
    cutoff_30 = (today - timedelta(days=30)).isoformat()
    cutoff_15 = (today - timedelta(days=15)).isoformat()
    now = datetime.utcnow().isoformat(timespec="seconds")

    log.info("Materializing cluster_summary (fast path); window cutoff = %s", cutoff_60)

    conn.execute("DELETE FROM cluster_summary")
    conn.commit()

    # ── Step 1: Bulk aggregate scalar fields ──────────────────────────────
    log.info("Step 1: bulk scalar aggregation...")
    conn.execute("""
        CREATE TEMP TABLE _cs AS
        SELECT
            pathogen, pds_acc,
            COUNT(*) AS n_total,
            SUM(CASE WHEN source_category='Human' THEN 1 ELSE 0 END) AS n_human,
            SUM(CASE WHEN source_category!='Human' THEN 1 ELSE 0 END) AS n_nonhuman,
            SUM(CASE WHEN source_category='Food' THEN 1 ELSE 0 END) AS n_food,
            SUM(CASE WHEN source_category='Animal' THEN 1 ELSE 0 END) AS n_animal,
            SUM(CASE WHEN source_category='Environment' THEN 1 ELSE 0 END) AS n_environment,
            MIN(collection_date) AS earliest_collection_date,
            MAX(collection_date) AS latest_collection_date,
            MIN(target_creation_date) AS earliest_target_creation_date,
            MAX(target_creation_date) AS latest_target_creation_date,
            SUM(CASE WHEN source_category='Human'
                     AND collection_date >= ? THEN 1 ELSE 0 END) AS new_humans_in_window,
            SUM(CASE WHEN source_category='Human'
                     AND collection_date >= ? THEN 1 ELSE 0 END) AS new_humans_30d,
            SUM(CASE WHEN source_category='Human'
                     AND collection_date >= ? THEN 1 ELSE 0 END) AS new_humans_15d
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
        GROUP BY pathogen, pds_acc
    """, (cutoff_60, cutoff_30, cutoff_15))
    conn.commit()

    total_clusters = conn.execute("SELECT COUNT(*) FROM _cs").fetchone()[0]
    log.info("Found %d distinct clusters", total_clusters)

    # ── Step 2: countries JSON per cluster ────────────────────────────────
    log.info("Step 2: country aggregation...")
    country_rows = conn.execute("""
        SELECT pds_acc, geo_country, COUNT(*) as n
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
          AND geo_country IS NOT NULL AND geo_country != ''
        GROUP BY pds_acc, geo_country
        ORDER BY pds_acc, n DESC
    """).fetchall()
    countries_by_pds: dict[str, list] = {}
    for r in country_rows:
        countries_by_pds.setdefault(r[0], []).append({"country": r[1], "n": r[2]})

    # ── Step 3: source summary JSON per cluster ───────────────────────────
    log.info("Step 3: source summary aggregation...")
    source_rows = conn.execute("""
        SELECT pds_acc, source_category,
               LOWER(TRIM(COALESCE(isolation_source,'(unspecified)'))) as src,
               COUNT(*) as n
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
        GROUP BY pds_acc, source_category, src
        ORDER BY pds_acc, n DESC
    """).fetchall()
    sources_by_pds: dict[str, list] = {}
    for r in source_rows:
        pds = r[0]
        if pds not in sources_by_pds:
            sources_by_pds[pds] = []
        if len(sources_by_pds[pds]) < 20:
            sources_by_pds[pds].append({"category": r[1], "source": r[2], "n": r[3]})

    # ── Step 4: host summary JSON per cluster ─────────────────────────────
    log.info("Step 4: host summary aggregation...")
    host_rows = conn.execute("""
        SELECT pds_acc, LOWER(TRIM(host)) as h, COUNT(*) as n
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
          AND source_category != 'Human'
          AND host IS NOT NULL AND host != ''
          AND LOWER(TRIM(host)) NOT IN ('missing','not collected','not provided','unknown','na')
        GROUP BY pds_acc, h
        ORDER BY pds_acc, n DESC
    """).fetchall()
    hosts_by_pds: dict[str, list] = {}
    for r in host_rows:
        pds = r[0]
        if pds not in hosts_by_pds:
            hosts_by_pds[pds] = []
        if len(hosts_by_pds[pds]) < 10:
            hosts_by_pds[pds].append({"host": r[1], "n": r[2]})

    # ── Step 5: recent human cases JSON per cluster ───────────────────────
    log.info("Step 5: recent human cases aggregation...")
    recent_human_rows = conn.execute("""
        SELECT pds_acc, pdt_acc, biosample_acc,
               collection_date, collection_date_raw,
               target_creation_date, geo_loc_name, geo_country
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
          AND source_category = 'Human'
          AND collection_date >= ?
        ORDER BY pds_acc, collection_date DESC, pdt_acc
    """, (cutoff_60,)).fetchall()
    recent_humans_by_pds: dict[str, list] = {}
    for r in recent_human_rows:
        recent_humans_by_pds.setdefault(r[0], []).append({
            "pdt": r[1], "biosample": r[2],
            "collection_date": r[4] or r[3],
            "date_added": r[5],
            "geo": r[6],
            "geo_country": r[7] if len(r) > 7 else "",
        })

    # ── Step 6: deposit lag per cluster ──────────────────────────────────
    log.info("Step 6: deposit lag aggregation...")
    lag_rows = conn.execute("""
        SELECT pds_acc,
               CAST(julianday(target_creation_date) - julianday(collection_date) AS INTEGER) as lag
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
          AND collection_date IS NOT NULL
          AND target_creation_date IS NOT NULL
          AND julianday(target_creation_date) >= julianday(collection_date)
        ORDER BY pds_acc, lag
    """).fetchall()
    lags_by_pds: dict[str, list] = {}
    for r in lag_rows:
        lags_by_pds.setdefault(r[0], []).append(r[1])

    # ── Step 7: histogram per cluster ─────────────────────────────────────
    log.info("Step 7: histogram aggregation...")
    hist_rows = conn.execute("""
        SELECT pds_acc,
               CAST(SUBSTR(collection_date,1,4) AS INTEGER) as yr,
               SUM(CASE WHEN source_category='Human' THEN 1 ELSE 0 END) as n_human,
               SUM(CASE WHEN source_category!='Human' THEN 1 ELSE 0 END) as n_nonhuman
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
          AND collection_date IS NOT NULL
          AND LENGTH(collection_date) >= 4
        GROUP BY pds_acc, yr
        ORDER BY pds_acc, yr
    """).fetchall()
    hist_by_pds: dict[str, list] = {}
    for r in hist_rows:
        hist_by_pds.setdefault(r[0], []).append(
            {"year": r[1], "n_human": r[2], "n_nonhuman": r[3]}
        )

    # ── Step 8: oldest isolates per cluster ───────────────────────────────
    log.info("Step 8: oldest isolate aggregation...")
    oldest_rows = conn.execute("""
        SELECT pds_acc, pdt_acc, biosample_acc, collection_date,
               collection_date_raw, geo_loc_name, geo_country, geo_admin1,
               isolation_source, source_category,
               ROW_NUMBER() OVER (
                   PARTITION BY pds_acc
                   ORDER BY collection_date ASC, pdt_acc ASC
               ) as rn_any,
               ROW_NUMBER() OVER (
                   PARTITION BY pds_acc, source_category
                   ORDER BY collection_date ASC, pdt_acc ASC
               ) as rn_cat
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
          AND collection_date IS NOT NULL
    """).fetchall()
    oldest_any_by_pds: dict = {}
    oldest_human_by_pds: dict = {}
    oldest_nonhuman_by_pds: dict = {}
    for r in oldest_rows:
        pds = r[0]
        pack = {
            "pdt": r[1], "biosample": r[2], "date": r[3],
            "date_raw": r[4], "geo": r[5], "geo_country": r[6],
            "geo_admin1": r[7], "source": r[8], "source_category": r[9],
        }
        if r[10] == 1:  # rn_any
            oldest_any_by_pds[pds] = pack
        if r[11] == 1:  # rn_cat
            if r[9] == "Human":
                oldest_human_by_pds[pds] = pack
            else:
                oldest_nonhuman_by_pds[pds] = pack

    # ── Step 9: latest assembly per cluster ───────────────────────────────
    log.info("Step 9: latest assembly aggregation...")
    asm_rows = conn.execute("""
        SELECT pds_acc, pdt_acc, biosample_acc, asm_acc,
               collection_date_raw, collection_date, target_creation_date
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
          AND asm_acc IS NOT NULL AND TRIM(asm_acc) != ''
        ORDER BY pds_acc, target_creation_date DESC, collection_date DESC, pdt_acc
    """).fetchall()
    latest_asm_by_pds: dict = {}
    for r in asm_rows:
        pds = r[0]
        if pds not in latest_asm_by_pds:
            latest_asm_by_pds[pds] = {
                "pdt": r[1], "biosample": r[2], "asm_acc": r[3],
                "collection_date": r[4] or r[5],
                "target_creation_date": r[6],
            }

    # ── Step 10: map locations per cluster ────────────────────────────────
    log.info("Step 10: map location aggregation...")
    map_rows = conn.execute("""
        SELECT pds_acc, geo_country, geo_admin1, geo_loc_name,
               source_category, collection_date, COUNT(*) as n
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
          AND collection_date IS NOT NULL
        GROUP BY pds_acc, geo_country, geo_admin1, source_category, collection_date
        ORDER BY pds_acc
    """).fetchall()
    map_members_by_pds: dict[str, list] = {}
    for r in map_rows:
        map_members_by_pds.setdefault(r[0], []).append({
            "geo_country": r[1], "geo_admin1": r[2], "geo_loc_name": r[3],
            "source_category": r[4], "collection_date": r[5], "n": r[6],
        })

    # ── Step 11: admin1 per cluster ───────────────────────────────────────
    log.info("Step 11: admin1 aggregation...")
    admin1_rows = conn.execute("""
        SELECT pds_acc, geo_country, geo_admin1, COUNT(*) as n
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
          AND geo_country IS NOT NULL
        GROUP BY pds_acc, geo_country, geo_admin1
    """).fetchall()
    admin1_by_pds: dict[str, dict] = {}
    for r in admin1_rows:
        pds, country, admin1, n = r[0], r[1], r[2], r[3]
        if pds not in admin1_by_pds:
            admin1_by_pds[pds] = {"by_country": {}, "unspecified": {}}
        if admin1:
            d = admin1_by_pds[pds]["by_country"].setdefault(country, {})
            d[admin1] = d.get(admin1, 0) + n
        else:
            admin1_by_pds[pds]["unspecified"][country] = (
                admin1_by_pds[pds]["unspecified"].get(country, 0) + n
            )

    # ── Step 12: signals + INSERT per cluster ─────────────────────────────
    log.info("Step 11b: bulk member fetch for signals computation...")
    members_by_pds: dict[str, list] = {}
    member_rows = conn.execute("""
        SELECT pdt_acc, pds_acc, pathogen, epi_type, source_category,
               geo_country, geo_admin1, geo_loc_name, isolation_source,
               collection_date, target_creation_date, food_origin,
               bioproject_acc, biosample_acc, host
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
    """).fetchall()
    for r in member_rows:
        members_by_pds.setdefault(r["pds_acc"], []).append(dict(r))
    log.info(f"  Fetched member data for {len(members_by_pds)} clusters")

    log.info("Step 12: computing signals and inserting rows...")
    from .. import signals as _signals
    from ..render import map as _map

    scalar_rows = conn.execute("""
        SELECT pathogen, pds_acc, n_total, n_human, n_nonhuman,
               n_food, n_animal, n_environment,
               earliest_collection_date, latest_collection_date,
               earliest_target_creation_date, latest_target_creation_date,
               new_humans_in_window, new_humans_30d, new_humans_15d
        FROM _cs
    """).fetchall()

    # Pre-fetch MLST STs
    st_rows = conn.execute(
        "SELECT pds_acc, mlst_st FROM cluster_typing"
    ).fetchall()
    mlst_st_by_pds = {r[0]: r[1] for r in st_rows}

    batch = []
    written = 0
    for sr in scalar_rows:
        pds = sr["pds_acc"]
        pathogen = sr["pathogen"]

        # Temporal span
        ec = sr["earliest_collection_date"]
        lc = sr["latest_collection_date"]
        temporal_span_days = None
        if ec and lc:
            try:
                temporal_span_days = (
                    date.fromisoformat(str(lc)[:10]) -
                    date.fromisoformat(str(ec)[:10])
                ).days
            except (ValueError, TypeError):
                pass

        # Deposit lag
        lags = sorted(lags_by_pds.get(pds, []))
        deposit_lag_median = lags[len(lags) // 2] if lags else None
        deposit_lag_mean = (sum(lags) // len(lags)) if lags else None

        # Histogram
        histogram = hist_by_pds.get(pds, [])
        histogram_max = max(
            (h["n_human"] + h["n_nonhuman"] for h in histogram), default=0
        )

        # Map locations
        map_locs = _map.aggregate_locations(
            map_members_by_pds.get(pds, []), today=today
        )

        # Admin1 cleanup
        raw_admin1 = admin1_by_pds.get(pds, {"by_country": {}, "unspecified": {}})
        admin1_data = {
            "by_country": {
                c: sorted(
                    [{"admin1": a, "n": n} for a, n in d.items()],
                    key=lambda x: -x["n"]
                )
                for c, d in raw_admin1["by_country"].items()
            },
            "unspecified": raw_admin1["unspecified"],
        }

        # Signals — use real member dicts
        signals_blob = _signals.compute_all_signals(
            members=members_by_pds.get(pds, []),
            amr_genes_per_isolate={},
            mlst_st=mlst_st_by_pds.get(pds),
            today=today,
            recent_window_days=window_days,
        )

        recent_humans = recent_humans_by_pds.get(pds, [])

        batch.append((
            pathogen, pds,
            sr["n_total"], sr["n_human"], sr["n_nonhuman"],
            sr["n_food"], sr["n_animal"], sr["n_environment"],
            ec, lc,
            sr["earliest_target_creation_date"],
            sr["latest_target_creation_date"],
            temporal_span_days,
            json.dumps(oldest_any_by_pds.get(pds)),
            json.dumps(oldest_human_by_pds.get(pds)),
            json.dumps(oldest_nonhuman_by_pds.get(pds)),
            json.dumps(histogram) if histogram else None,
            histogram_max,
            deposit_lag_median, deposit_lag_mean,
            json.dumps(hosts_by_pds.get(pds)) if hosts_by_pds.get(pds) else None,
            json.dumps(map_locs) if map_locs else None,
            json.dumps(latest_asm_by_pds.get(pds)) if latest_asm_by_pds.get(pds) else None,
            json.dumps(signals_blob),
            json.dumps(countries_by_pds.get(pds, [])),
            json.dumps(admin1_data),
            json.dumps(sources_by_pds.get(pds, [])),
            sr["new_humans_in_window"],
            sr["new_humans_30d"],
            sr["new_humans_15d"],
            json.dumps([r["pdt"] for r in recent_humans]),
            json.dumps(recent_humans),
            window_days,
            now,
        ))
        written += 1

        if len(batch) >= 500:
            conn.executemany("""
                INSERT INTO cluster_summary (
                    pathogen, pds_acc,
                    n_total, n_human, n_nonhuman, n_food, n_animal, n_environment,
                    earliest_collection_date, latest_collection_date,
                    earliest_target_creation_date, latest_target_creation_date,
                    temporal_span_days,
                    oldest_isolate_json, oldest_human_json, oldest_nonhuman_json,
                    histogram_json, histogram_max_year_count,
                    deposit_lag_median, deposit_lag_mean,
                    host_summary_json,
                    map_locations_json, latest_assembly_json,
                    signals_json,
                    countries_json, admin1_json, source_summary_json,
                    new_humans_in_window, new_humans_30d, new_humans_15d,
                    new_humans_in_window_pdts_json,
                    new_humans_in_window_dates_json,
                    window_days, refreshed_at
                ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """, batch)
            conn.commit()
            log.info("  Inserted %d / %d clusters...", written, total_clusters)
            batch = []

    if batch:
        conn.executemany("""
            INSERT INTO cluster_summary (
                pathogen, pds_acc,
                n_total, n_human, n_nonhuman, n_food, n_animal, n_environment,
                earliest_collection_date, latest_collection_date,
                earliest_target_creation_date, latest_target_creation_date,
                temporal_span_days,
                oldest_isolate_json, oldest_human_json, oldest_nonhuman_json,
                histogram_json, histogram_max_year_count,
                deposit_lag_median, deposit_lag_mean,
                host_summary_json,
                map_locations_json, latest_assembly_json,
                signals_json,
                countries_json, admin1_json, source_summary_json,
                new_humans_in_window, new_humans_30d, new_humans_15d,
                new_humans_in_window_pdts_json,
                new_humans_in_window_dates_json,
                window_days, refreshed_at
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, batch)
        conn.commit()

    conn.execute("DROP TABLE IF EXISTS _cs")
    conn.commit()

    log.info("Materialized %d cluster_summary rows", written)
    with_humans = conn.execute(
        "SELECT COUNT(*) FROM cluster_summary WHERE n_human > 0"
    ).fetchone()[0]
    mixed = conn.execute(
        "SELECT COUNT(*) FROM cluster_summary WHERE n_human > 0 AND n_nonhuman > 0"
    ).fetchone()[0]
    recent = conn.execute(
        "SELECT COUNT(*) FROM cluster_summary WHERE new_humans_in_window > 0"
    ).fetchone()[0]
    log.info(
        "Cluster shape: %d with humans; %d mixed; %d with recent human activity (last %dd)",
        with_humans, mixed, recent, window_days,
    )
    return written
