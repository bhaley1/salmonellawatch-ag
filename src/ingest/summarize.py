"""Materialize the cluster_summary table.

After every ingest, rebuild the per-cluster summary that the renderer uses.
This is the central join that makes the dashboard query fast: instead of
recomputing aggregates on every page load, we precompute them once per run.

Specifically: for each (pathogen, pds_acc), we compute:
  - Total counts: total, human, nonhuman, food, animal, environment
  - Date ranges: earliest/latest collection and target_creation
  - Geography: countries seen, sorted by frequency
  - Source summary: top nonhuman sources
  - Recent activity: new human PDTs in the last 30 days, with their dates/geo
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
    """Rebuild cluster_summary table. Returns number of clusters written."""
    window_days = window_days or config.RECENT_WINDOW_DAYS
    today = today or date.today()
    cutoff = (today - timedelta(days=window_days)).isoformat()
    now = datetime.utcnow().isoformat(timespec="seconds")

    log.info("Materializing cluster_summary; window cutoff = %s", cutoff)

    # Clear the table first — it's small, fast to rebuild from scratch
    conn.execute("DELETE FROM cluster_summary")
    conn.commit()

    # Find all (pathogen, pds_acc) combos that have at least one isolate.
    # Drop NULL pds_acc (isolates not in any cluster — not surveillance signal).
    clusters = conn.execute("""
        SELECT pathogen, pds_acc, COUNT(*) AS n_total
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
        GROUP BY pathogen, pds_acc
    """).fetchall()

    log.info("Found %d distinct clusters", len(clusters))

    written = 0
    for cluster_row in clusters:
        pathogen = cluster_row["pathogen"]
        pds_acc = cluster_row["pds_acc"]

        # Pull all members of this cluster
        members = conn.execute("""
            SELECT pdt_acc, epi_type, source_category, geo_country, geo_admin1,
                   geo_loc_name, isolation_source, collection_date,
                   collection_date_raw, target_creation_date,
                   serovar, asm_acc, biosample_acc, food_origin, ifsac_category,
                   host_disease, host, bioproject_acc
            FROM isolates
            WHERE pathogen = ? AND pds_acc = ?
        """, (pathogen, pds_acc)).fetchall()

        # Aggregate counts
        n_human = sum(1 for m in members if m["source_category"] == "Human")
        n_food = sum(1 for m in members if m["source_category"] == "Food")
        n_animal = sum(1 for m in members if m["source_category"] == "Animal")
        n_environment = sum(1 for m in members if m["source_category"] == "Environment")
        n_nonhuman = len(members) - n_human

        # Date ranges
        coll_dates = [m["collection_date"] for m in members if m["collection_date"]]
        tgt_dates = [m["target_creation_date"] for m in members if m["target_creation_date"]]
        earliest_coll = min(coll_dates) if coll_dates else None
        latest_coll = max(coll_dates) if coll_dates else None
        earliest_tgt = min(tgt_dates) if tgt_dates else None
        latest_tgt = max(tgt_dates) if tgt_dates else None

        # Temporal span: how long the cluster has been "active" by collection date.
        # A short span (weeks-months) suggests a discrete outbreak; a long span
        # (years) suggests persistent environmental contamination. For Listeria
        # specifically, processing-facility persistence over 5-20 years is well
        # documented, so this is a critical interpretive signal.
        temporal_span_days: int | None = None
        if earliest_coll and latest_coll:
            from datetime import date as _date
            try:
                e = _date.fromisoformat(earliest_coll) if isinstance(earliest_coll, str) else earliest_coll
                l = _date.fromisoformat(latest_coll) if isinstance(latest_coll, str) else latest_coll
                temporal_span_days = (l - e).days
            except (ValueError, TypeError):
                temporal_span_days = None

        # Oldest-isolate signatures — three separate ones because the
        # epidemiological meaning of each is different:
        #   oldest of any kind = when did the signal first appear at all
        #   oldest human       = when did the cluster first cause disease
        #   oldest nonhuman    = when did the (potential) source first appear
        # Among ties, prefer rows with more specific geography (admin1 over country).
        def _oldest_isolate_pack(m) -> dict:
            return {
                "pdt": m["pdt_acc"],
                "biosample": m["biosample_acc"],
                "date": m["collection_date"],
                "date_raw": m["collection_date_raw"],
                "geo": m["geo_loc_name"],
                "geo_country": m["geo_country"],
                "geo_admin1": m["geo_admin1"],
                "source": m["isolation_source"],
                "source_category": m["source_category"],
            }

        def _earliest(rows) -> dict | None:
            dated = [r for r in rows if r["collection_date"]]
            if not dated:
                return None
            # Sort by (date asc, admin1 specificity desc as tiebreaker, pdt for stability)
            dated.sort(key=lambda r: (
                r["collection_date"],
                0 if r["geo_admin1"] else 1,
                r["pdt_acc"],
            ))
            return _oldest_isolate_pack(dated[0])

        oldest_any = _earliest(list(members))
        oldest_human = _earliest([m for m in members if m["source_category"] == "Human"])
        oldest_nonhuman = _earliest([m for m in members if m["source_category"] != "Human"])

        oldest_isolate_json = json.dumps(oldest_any) if oldest_any else None
        oldest_human_json = json.dumps(oldest_human) if oldest_human else None
        oldest_nonhuman_json = json.dumps(oldest_nonhuman) if oldest_nonhuman else None

        # Collection-to-deposit lag — the typical delay between when an
        # isolate is collected from the world and when it appears in NCBI.
        # Surveillance professionals use this to interpret the freshness of
        # the dashboard data: a 7-day lag means we're seeing near-real-time
        # surveillance; a 180-day lag means most "new" deposits reflect old
        # cases being uploaded. Computed only on members with both dates.
        lag_values: list[int] = []
        from datetime import date as _date_cls
        for m in members:
            cd = m["collection_date"]
            tcd = m["target_creation_date"]
            if not cd or not tcd:
                continue
            try:
                cdd = _date_cls.fromisoformat(str(cd)[:10])
                tcdd = _date_cls.fromisoformat(str(tcd)[:10])
                lag_days = (tcdd - cdd).days
                if lag_days >= 0:  # negative lags are submission-date metadata errors
                    lag_values.append(lag_days)
            except (ValueError, TypeError):
                continue
        deposit_lag_median: int | None = None
        deposit_lag_mean: int | None = None
        if lag_values:
            sorted_lags = sorted(lag_values)
            deposit_lag_median = sorted_lags[len(sorted_lags) // 2]
            deposit_lag_mean = sum(lag_values) // len(lag_values)

        # Histogram — per-year human vs. nonhuman counts.
        # Only counts isolates with parseable collection_date. Year is the
        # first 4 chars of the ISO date string.
        year_counts: dict[int, dict[str, int]] = {}
        for m in members:
            cd = m["collection_date"]
            if not cd:
                continue
            try:
                year = int(str(cd)[:4])
            except (ValueError, TypeError):
                continue
            yc = year_counts.setdefault(year, {"human": 0, "nonhuman": 0})
            if m["source_category"] == "Human":
                yc["human"] += 1
            else:
                yc["nonhuman"] += 1
        # Sort years ascending and serialize as a stable list
        histogram = [
            {"year": yr, "n_human": year_counts[yr]["human"], "n_nonhuman": year_counts[yr]["nonhuman"]}
            for yr in sorted(year_counts.keys())
        ]
        histogram_json = json.dumps(histogram) if histogram else None
        histogram_max_year_count = max(
            (h["n_human"] + h["n_nonhuman"] for h in histogram), default=0,
        ) if histogram else 0

        # Map locations — full cluster geographic footprint.
        #
        # The map shows ALL isolates with a known collection_date and a
        # known location, regardless of how old they are. For long-persistent
        # Listeria clusters (many span 10-30+ years) a 1-year window misses
        # the cluster's actual geographic distribution. We do encode recency
        # via dot opacity at render time: recent dots full-opacity, old dots
        # progressively faded.
        map_members = [
            dict(m) for m in members
            if m["collection_date"]
        ]
        from ..render import map as _map  # local import to avoid circular
        map_locations = _map.aggregate_locations(map_members, today=today)
        map_locations_json = json.dumps(map_locations) if map_locations else None

        # Latest-assembled isolate: the most-recently-deposited isolate in
        # this cluster that has an asm_acc. Used by the dashboard to give
        # users a "view the latest assembled genome for this cluster" link.
        #
        # Sort by target_creation_date DESC (when NCBI processed it), then
        # collection_date DESC as tiebreak. Only consider isolates with
        # a populated asm_acc — many isolates are SNP-typed without going
        # through assembly, so ~50% of NCBI PD isolates have asm_acc=None.
        assembled = [
            m for m in members
            if m["asm_acc"] and m["asm_acc"].strip()
        ]
        latest_assembly_json: str | None = None
        if assembled:
            assembled.sort(
                key=lambda m: (
                    m["target_creation_date"] or "",
                    m["collection_date"] or "",
                    m["pdt_acc"],
                ),
                reverse=True,
            )
            la = assembled[0]
            latest_assembly_json = json.dumps({
                "pdt": la["pdt_acc"],
                "biosample": la["biosample_acc"],
                "asm_acc": la["asm_acc"],
                "collection_date": la["collection_date_raw"] or la["collection_date"],
                "target_creation_date": la["target_creation_date"],
            })

        # Countries, sorted by count desc
        country_counts: dict[str, int] = {}
        for m in members:
            c = m["geo_country"]
            if c:
                country_counts[c] = country_counts.get(c, 0) + 1
        countries = sorted(country_counts.items(), key=lambda x: -x[1])
        countries_json = json.dumps([{"country": c, "n": n} for c, n in countries])

        # Admin1 breakdown (state/region within country), grouped by country.
        # Stored as {country: [{admin1: str, n: int}, ...], "_unspecified": {country: n}}
        # The renderer uses this to show specific geography (e.g. "USA — Maryland (3),
        # New York (2)") instead of bare country counts.
        admin1_by_country: dict[str, dict[str, int]] = {}
        unspecified_by_country: dict[str, int] = {}
        for m in members:
            country = m["geo_country"]
            admin1 = m["geo_admin1"]
            if not country:
                continue
            if admin1:
                d = admin1_by_country.setdefault(country, {})
                d[admin1] = d.get(admin1, 0) + 1
            else:
                unspecified_by_country[country] = unspecified_by_country.get(country, 0) + 1
        admin1_data = {
            "by_country": {
                c: sorted([{"admin1": a, "n": n} for a, n in d.items()],
                          key=lambda x: -x["n"])
                for c, d in admin1_by_country.items()
            },
            "unspecified": unspecified_by_country,
        }
        admin1_json = json.dumps(admin1_data)

        # Source summary — top sources across ALL categories (human + nonhuman).
        # This was previously nonhuman-only. We now include human-source strings
        # (e.g. 'blood', 'CSF') because they're epidemiologically meaningful
        # alongside food/animal/environment counts. Renderer can filter if needed.
        # 'Unknown' source category is included to make the totals honest.
        source_counts: dict[tuple[str, str], int] = {}
        for m in members:
            cat = m["source_category"] or "Unknown"
            src = (m["isolation_source"] or "").strip().lower() or "(unspecified)"
            key = (cat, src)
            source_counts[key] = source_counts.get(key, 0) + 1
        sources_sorted = sorted(source_counts.items(), key=lambda x: -x[1])[:20]
        source_summary_json = json.dumps([
            {"category": cat, "source": src, "n": n}
            for (cat, src), n in sources_sorted
        ])

        # Host species summary — actual host strings from animal isolates.
        # Skips humans (which are already collapsed into 'Human' source category)
        # and skips blank/unknown hosts.
        #
        # NCBI host strings are inconsistent: "Bos taurus", "bovine", "cattle",
        # "cow" all mean the same thing. We don't try to canonicalize because
        # the canonical mapping is contested in veterinary literature. We just
        # display the strings as submitted and let the reader collapse.
        host_counts: dict[str, int] = {}
        for m in members:
            sc = m["source_category"]
            if sc == "Human" or sc == "Unknown":
                continue
            h = (m["host"] or "").strip().lower()
            if not h or h in ("missing", "not collected", "not provided", "unknown", "na"):
                continue
            host_counts[h] = host_counts.get(h, 0) + 1
        hosts_sorted = sorted(host_counts.items(), key=lambda x: -x[1])[:10]
        host_summary_json = json.dumps([
            {"host": h, "n": n} for h, n in hosts_sorted
        ]) if hosts_sorted else None

        # Recent human cases (the surveillance signal).
        #
        # Filter on collection_date, NOT target_creation_date. This means
        # "actually recent illness" rather than "isolates NCBI happened to
        # ingest recently regardless of when sampled." Backlogged deposits
        # — labs uploading 2019 sequences in 2026 — are valid surveillance
        # data but they're not new cases and shouldn't appear here.
        #
        # We also compute counts at 15d and 30d windows for the header
        # recency gradient ("9 in 60d · 4 in 30d · 1 in 15d").
        cutoff_60 = cutoff   # already computed above using window_days
        cutoff_30 = (today - timedelta(days=30)).isoformat()
        cutoff_15 = (today - timedelta(days=15)).isoformat()

        def _has_recent_collection(m, c):
            cd = m["collection_date"]
            return bool(cd) and str(cd) >= c

        recent_humans = [
            m for m in members
            if m["source_category"] == "Human"
            and _has_recent_collection(m, cutoff_60)
        ]
        # Sort by collection_date desc, pdt_acc tiebreak
        recent_humans.sort(
            key=lambda m: (m["collection_date"] or "", m["pdt_acc"]),
            reverse=True,
        )
        new_humans_n = len(recent_humans)
        new_humans_30d = sum(
            1 for m in members
            if m["source_category"] == "Human" and _has_recent_collection(m, cutoff_30)
        )
        new_humans_15d = sum(
            1 for m in members
            if m["source_category"] == "Human" and _has_recent_collection(m, cutoff_15)
        )

        new_humans_pdts = [m["pdt_acc"] for m in recent_humans]
        new_humans_pdts_json = json.dumps(new_humans_pdts)
        new_humans_dates_json = json.dumps([
            {
                "pdt": m["pdt_acc"],
                "biosample": m["biosample_acc"],
                # Keep date_added in the data export for transparency, but
                # the template no longer displays it.
                "date_added": m["target_creation_date"],
                "collection_date": m["collection_date_raw"] or m["collection_date"],
                "geo": m["geo_loc_name"],
            }
            for m in recent_humans
        ])

        # Per-cluster derived signals (geographic spread, import, AMR, etc.)
        # AMR data must be loaded into a {pdt_acc: {gene_symbol, ...}} dict.
        amr_rows = conn.execute("""
            SELECT a.pdt_acc, a.gene_symbol
            FROM isolate_amr a
            JOIN isolates i ON a.pdt_acc = i.pdt_acc
            WHERE i.pathogen = ? AND i.pds_acc = ?
        """, (pathogen, pds_acc)).fetchall()
        amr_by_pdt: dict[str, set[str]] = {}
        for r in amr_rows:
            amr_by_pdt.setdefault(r["pdt_acc"], set()).add(r["gene_symbol"])

        # MLST ST for this cluster (joined from cluster_typing if present)
        ct_row = conn.execute(
            "SELECT mlst_st FROM cluster_typing WHERE pathogen = ? AND pds_acc = ?",
            (pathogen, pds_acc),
        ).fetchone()
        mlst_st = ct_row["mlst_st"] if ct_row else None

        # Convert sqlite.Row members to plain dicts for signals.py
        member_dicts = [dict(m) for m in members]

        from .. import signals as _signals
        signals_blob = _signals.compute_all_signals(
            members=member_dicts,
            amr_genes_per_isolate=amr_by_pdt,
            mlst_st=mlst_st,
            today=today,
            recent_window_days=window_days,
        )
        signals_json = json.dumps(signals_blob)

        conn.execute("""
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
                window_days,
                refreshed_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            pathogen, pds_acc,
            len(members), n_human, n_nonhuman, n_food, n_animal, n_environment,
            earliest_coll, latest_coll, earliest_tgt, latest_tgt,
            temporal_span_days,
            oldest_isolate_json, oldest_human_json, oldest_nonhuman_json,
            histogram_json, histogram_max_year_count,
            deposit_lag_median, deposit_lag_mean,
            host_summary_json,
            map_locations_json, latest_assembly_json,
            signals_json,
            countries_json, admin1_json, source_summary_json,
            new_humans_n, new_humans_30d, new_humans_15d,
            new_humans_pdts_json, new_humans_dates_json,
            window_days,
            now,
        ))
        written += 1

    conn.commit()
    log.info("Materialized %d cluster_summary rows", written)

    # Bonus diagnostics for the log
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
