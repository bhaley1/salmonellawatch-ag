"""Phase 1 minimal renderer.

Just enough to prove the SQLite → HTML pipeline works. The full UX from
the design doc (progressive disclosure, plain language, tooltips, beta
banner) ships in Phase 3.

For now we render exactly one page (the recent activity dashboard) plus
a JSON export. Both pull from cluster_summary via the queries module.
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from .. import config, db
from . import queries

log = logging.getLogger("pathogen-watch-render")


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def render() -> int:
    site_dir = config.SITE_DIR
    site_dir.mkdir(parents=True, exist_ok=True)

    env = Environment(
        loader=FileSystemLoader(str(config.TEMPLATES_DIR)),
        autoescape=select_autoescape(["html", "xml"]),
        trim_blocks=True,
        lstrip_blocks=True,
    )
    # Register prose helpers as template filters
    from . import histogram, map as cluster_map, prose
    env.filters["span_label"] = prose.format_temporal_span
    env.filters["span_interp"] = prose.span_interpretation
    env.filters["geography_summary"] = prose.format_geography
    env.filters["oldest_isolate_label"] = prose.format_oldest_isolate
    env.globals["ncbi_link_for"] = prose.ncbi_link_for

    def _hist(cluster):
        return histogram.render_histogram_svg(
            cluster.get("histogram") or [],
            cluster.get("histogram_max_year_count"),
        )
    env.filters["histogram_svg"] = _hist

    def _map(cluster):
        return cluster_map.render_cluster_map_svg(cluster.get("map_locations") or [])
    env.filters["map_svg"] = _map

    with db.connect() as conn:
        clusters, serotypes_for_dropdown = queries.get_recent_activity_clusters(conn)
        totals = queries.get_totals(conn)
        pathogen_counts = queries.get_pathogen_counts(conn)
        latest_release = queries.get_latest_release(conn)

    log.info(
        "Rendering: %d active clusters, %d total clusters, latest=%s",
        len(clusters),
        totals.get("n_clusters", 0),
        (latest_release or {}).get("pdg_release"),
    )

    # Compute the set of distinct countries that have at least one human
    # isolate in some active cluster. These become the dropdown filter options.
    # Sort with USA first (most users care about US clusters), then
    # alphabetical for the rest. This is the data-driven approach: countries
    # not present in active clusters don't appear in the dropdown.
    country_set: set[str] = set()
    for c in clusters:
        signals = c.get("signals") or {}
        imp = signals.get("import_signal") or {}
        for country in imp.get("human_countries") or []:
            country_set.add(country)
    # Build the sorted dropdown list
    countries_for_dropdown = []
    if "United States" in country_set:
        countries_for_dropdown.append("United States")
        country_set.remove("United States")
    countries_for_dropdown.extend(sorted(country_set))

    # Annotate each cluster with its human countries as a comma-separated
    # attribute string for the data-human-countries HTML attribute.
    for c in clusters:
        signals = c.get("signals") or {}
        imp = signals.get("import_signal") or {}
        c["human_countries_attr"] = ",".join(imp.get("human_countries") or [])

    common = {
        "generated_at": datetime.utcnow().strftime("%B %-d, %Y at %H:%M UTC"),
        "totals": totals,
        "pathogen_counts": pathogen_counts,
        "latest_release": latest_release,
        "window_days": config.RECENT_WINDOW_DAYS,
        "countries_for_dropdown": countries_for_dropdown,
        "review_mode": config.REVIEW_MODE,
    }

    # ---- Render index.html ----
    soc_serotypes = {
        'Enteritidis', 'Typhimurium', 'I 4,[5],12:i:-', 'Heidelberg',
        'Infantis', 'Newport', 'Uganda', 'Braenderup', 'Muenchen',
        'Montevideo', 'Javiana', 'Reading', 'Dublin', 'Oranienburg',
        'Potsdam', 'Thompson', 'Saintpaul', 'Hadar', 'Schwarzengrund',
        'Anatum', 'Berta',
    }
    rendered = env.get_template("index.html").render(
        clusters=clusters,
        serotypes_for_dropdown=serotypes_for_dropdown,
        soc_serotypes=soc_serotypes,
        **common,
    )
    (site_dir / "index.html").write_text(rendered)
    log.info("Wrote index.html (%d clusters)", len(clusters))

    # ---- Render methods.html ----
    methods_rendered = env.get_template("methods.html").render(**common)
    (site_dir / "methods.html").write_text(methods_rendered)
    log.info("Wrote methods.html")

    # ---- Data exports ----
    data_dir = site_dir / "data"
    data_dir.mkdir(exist_ok=True)
    (data_dir / "active-clusters.json").write_text(
        json.dumps(clusters, indent=2, default=str)
    )
    log.info("Wrote active-clusters.json")

    # ---- Static assets ----
    assets_src = Path(__file__).parent.parent.parent / "templates" / "assets"
    assets_dst = site_dir / "assets"
    if assets_src.exists():
        if assets_dst.exists():
            shutil.rmtree(assets_dst)
        shutil.copytree(assets_src, assets_dst)
        log.info("Copied static assets")

    # ---- GitHub Pages: .nojekyll ----
    (site_dir / ".nojekyll").write_text("")

    log.info("Render complete -> %s", site_dir)
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Pathogen Watch v2 render")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    setup_logging(args.verbose)
    return render()


if __name__ == "__main__":
    sys.exit(main())
