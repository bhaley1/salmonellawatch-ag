"""Inline SVG world map with split-dot isolate markers.

Renders a small per-cluster world map showing where isolates have been
collected in a 1-year window. Each dot is a split circle:
  - red half = human isolate portion at that location
  - teal half = nonhuman isolate portion at that location
  - dot radius scales with sqrt(total count)
  - SVG <title> elements provide hover tooltips with full breakdowns

Locations aggregate by:
  - US state when geo_country is United States and geo_admin1 is a US state
  - Country otherwise

No JavaScript. Server-rendered. ~5-15KB per cluster card.

Projection: simple equirectangular (lon→x, lat→y). Not strictly correct
for areas, but for a small inset showing dots, perfectly adequate.
"""

from __future__ import annotations

import math
from typing import Iterable

from ..lookups import centroids


# Canvas dimensions. Designed to fit in the right column of a cluster card.
W = 320
H = 180

# Padding around the projected area
PAD_X = 6
PAD_Y = 6

# Map extent: longitude -180..180, latitude -60..80 (cuts off Antarctica, much
# of the Arctic to better use the canvas)
LON_MIN, LON_MAX = -180.0, 180.0
LAT_MIN, LAT_MAX = -60.0, 80.0

# Visual constants
COLOR_HUMAN = "#b91c1c"        # rust — matches the rest of the dashboard
COLOR_NONHUMAN = "#0f766e"     # teal
COLOR_OUTLINE = "#fff"         # white outline around dot for contrast
COLOR_GRATICULE = "#e2e8f0"    # faint gridlines, visible over both land and ocean
COLOR_LAND = "#ebe6d4"         # warm off-white, slightly darker than page
COLOR_LAND_STROKE = "#cdc6ad"  # darker edge so land masses read clearly
COLOR_OCEAN = "#f5f3ec"        # very subtle blue-tinted off-white

# Dot size scaling: radius = MIN_R + sqrt(count) * SCALE, clamped to MAX_R.
MIN_R = 2.5
MAX_R = 11.0
SCALE = 1.6


def _project(lon: float, lat: float) -> tuple[float, float]:
    """Equirectangular projection. Returns (x, y) in SVG coords."""
    # x = (lon - LON_MIN) / (LON_MAX - LON_MIN) * (W - 2 PAD_X) + PAD_X
    x = PAD_X + (lon - LON_MIN) / (LON_MAX - LON_MIN) * (W - 2 * PAD_X)
    # y = (LAT_MAX - lat) / (LAT_MAX - LAT_MIN) * (H - 2 PAD_Y) + PAD_Y
    y = PAD_Y + (LAT_MAX - lat) / (LAT_MAX - LAT_MIN) * (H - 2 * PAD_Y)
    return (x, y)


def _radius_for(total_count: int) -> float:
    """Scale dot radius from isolate count using sqrt-scaling (area ≈ count)."""
    if total_count <= 0:
        return MIN_R
    r = MIN_R + math.sqrt(total_count) * SCALE
    return min(MAX_R, r)


def _opacity_for(days_ago: int | None) -> float:
    """Fade dots based on how old their most-recent isolate is.

    Recent (≤365 days):  full opacity 1.0
    Old (≥5 years):      faded to 0.35 (still visible but de-emphasized)
    Between:             linear interpolation
    None (unknown date): 0.55 (medium fade)
    """
    if days_ago is None:
        return 0.55
    if days_ago <= 365:
        return 1.0
    if days_ago >= 1825:  # 5 years
        return 0.35
    # Linear interp from 1.0 at 365 days to 0.35 at 1825 days
    frac = (days_ago - 365) / (1825 - 365)
    return 1.0 - frac * (1.0 - 0.35)


def _split_dot_svg(
    cx: float, cy: float, r: float,
    n_human: int, n_nonhuman: int,
    location_label: str,
    opacity: float = 1.0,
) -> str:
    """One split-dot. The left half is human, right half is nonhuman.

    When the location has only one source type, render a solid circle of
    that color (cleaner than a half-circle filling).

    Implementation: render the two halves as SVG paths so the proportions
    don't have to be 50/50. The vertical split position is determined by
    the ratio of human/(human+nonhuman).
    """
    total = n_human + n_nonhuman
    if total == 0:
        return ""

    title_text = f"{location_label}: {n_human} human, {n_nonhuman} nonhuman"
    opa = f' opacity="{opacity:.2f}"' if opacity < 1.0 else ""

    if n_human == 0:
        # Pure nonhuman — solid teal circle
        return (
            f'<g class="map-dot map-dot-nonhuman"{opa}>'
            f'<title>{title_text}</title>'
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
            f'fill="{COLOR_NONHUMAN}" stroke="{COLOR_OUTLINE}" stroke-width="0.8"/>'
            f'</g>'
        )
    if n_nonhuman == 0:
        # Pure human — solid red circle
        return (
            f'<g class="map-dot map-dot-human"{opa}>'
            f'<title>{title_text}</title>'
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
            f'fill="{COLOR_HUMAN}" stroke="{COLOR_OUTLINE}" stroke-width="0.8"/>'
            f'</g>'
        )

    # Mixed — vertical split with two paths
    left_path = (
        f'M {cx:.1f} {cy - r:.1f} '
        f'A {r:.1f} {r:.1f} 0 0 0 {cx:.1f} {cy + r:.1f} Z'
    )
    right_path = (
        f'M {cx:.1f} {cy - r:.1f} '
        f'A {r:.1f} {r:.1f} 0 0 1 {cx:.1f} {cy + r:.1f} Z'
    )

    return (
        f'<g class="map-dot map-dot-mixed"{opa}>'
        f'<title>{title_text}</title>'
        f'<path d="{left_path}" fill="{COLOR_HUMAN}"/>'
        f'<path d="{right_path}" fill="{COLOR_NONHUMAN}"/>'
        f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" '
        f'fill="none" stroke="{COLOR_OUTLINE}" stroke-width="0.8"/>'
        f'</g>'
    )


def _continents_svg() -> str:
    """Filled continent silhouettes as the map's land background.

    Each continent ring is projected and rendered as an SVG <path> with
    a filled, slightly off-white fill so it reads as land vs. ocean.
    Drawn once per map (cheap — total polygon count is ~80 vertices).
    """
    from . import continents
    parts = []
    for name, poly in continents.all_polygons():
        if not poly:
            continue
        # Project the ring
        proj = [_project(lon, lat) for lon, lat in poly]
        if not proj:
            continue
        # Build the SVG path: M first, L for each subsequent vertex, Z to close
        first = proj[0]
        d_parts = [f"M{first[0]:.1f},{first[1]:.1f}"]
        for x, y in proj[1:]:
            d_parts.append(f"L{x:.1f},{y:.1f}")
        d_parts.append("Z")
        parts.append(
            f'<path d="{"".join(d_parts)}" '
            f'fill="{COLOR_LAND}" stroke="{COLOR_LAND_STROKE}" '
            f'stroke-width="0.4" stroke-linejoin="round"/>'
        )
    return "".join(parts)


def _graticule_svg() -> str:
    """Faint graticule (lon/lat gridlines).

    Drawn after continents so the gridlines sit lightly on top of land.
    Lighter than the dots so they don't compete visually.
    """
    lines = []
    # Vertical lines: every 30 degrees of longitude
    for lon in range(int(LON_MIN), int(LON_MAX) + 1, 30):
        x, _ = _project(lon, 0)
        lines.append(
            f'<line x1="{x:.1f}" y1="{PAD_Y}" x2="{x:.1f}" y2="{H - PAD_Y}" '
            f'stroke="{COLOR_GRATICULE}" stroke-width="0.4" opacity="0.6"/>'
        )
    # Horizontal lines: every 30 degrees of latitude
    for lat in range(int(LAT_MIN), int(LAT_MAX) + 1, 30):
        _, y = _project(0, lat)
        lines.append(
            f'<line x1="{PAD_X}" y1="{y:.1f}" x2="{W - PAD_X}" y2="{y:.1f}" '
            f'stroke="{COLOR_GRATICULE}" stroke-width="0.4" opacity="0.6"/>'
        )
    return "".join(lines)


def _ocean_svg() -> str:
    """Background rectangle representing ocean."""
    return f'<rect x="0" y="0" width="{W}" height="{H}" fill="{COLOR_OCEAN}"/>'


def render_cluster_map_svg(locations: list[dict]) -> str:
    """Render a cluster's geographic footprint as a small SVG world map.

    `locations` is a list of dicts:
        {
          "label": str,              # display name (e.g. "USA - Maryland")
          "lon": float, "lat": float,
          "n_human": int, "n_nonhuman": int,
        }

    Dots are drawn smallest-first so large dots don't get hidden behind
    small adjacent ones — common in dense regions like Western Europe.
    """
    if not locations:
        return _empty_svg()

    background = _ocean_svg() + _continents_svg() + _graticule_svg()

    # Sort by total count ASC so big dots render on top
    sorted_locs = sorted(
        locations,
        key=lambda l: (l["n_human"] + l["n_nonhuman"]),
    )

    dots = []
    for loc in sorted_locs:
        x, y = _project(loc["lon"], loc["lat"])
        r = _radius_for(loc["n_human"] + loc["n_nonhuman"])
        opacity = _opacity_for(loc.get("days_ago"))
        dots.append(_split_dot_svg(
            x, y, r,
            loc["n_human"], loc["n_nonhuman"],
            loc["label"],
            opacity=opacity,
        ))

    svg = (
        f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" '
        f'class="cluster-map" role="img" '
        f'aria-label="World map showing isolate locations: red dots are human cases, '
        f'teal dots are nonhuman sources, split dots are mixed locations">'
        f'{background}'
        f'{"".join(dots)}'
        f'</svg>'
    )
    return svg


def _empty_svg() -> str:
    background = _ocean_svg() + _continents_svg() + _graticule_svg()
    return (
        f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" '
        f'class="cluster-map cluster-map-empty">'
        f'{background}'
        f'<text x="{W/2}" y="{H/2}" font-size="11" fill="#777" '
        f'text-anchor="middle" dominant-baseline="middle" font-style="italic">'
        f'No mapped isolates'
        f'</text></svg>'
    )


def aggregate_locations(
    members: list[dict],
    today: "date | None" = None,
) -> list[dict]:
    """Aggregate cluster member dicts into per-location split-dot data.

    Rules:
      - US isolates with a known state → aggregate by state (state centroid)
      - All other isolates → aggregate by country (country centroid)
      - Isolates without a geocodable location are dropped silently
      - Members must be pre-filtered to the time window the caller wants

    Each location bucket carries a most_recent_days_ago field so the
    renderer can fade old dots while keeping recent ones at full opacity.

    Returns the list of location dicts the renderer wants.
    """
    # Local imports to avoid circular dependencies
    from datetime import date as _date
    from ..lookups import geography as geo_lookup

    if today is None:
        today = _date.today()

    # Keys: (canonical_country, admin1_or_None_for_country) → counts
    by_key: dict[tuple[str, str | None], dict] = {}

    for m in members:
        country_raw = (m.get("geo_country") or "").strip()
        country = geo_lookup.canonical_country(country_raw) if country_raw else None
        admin1 = (m.get("geo_admin1") or "").strip() or None
        if not country:
            continue

        # Resolve centroid + label
        if country == "United States" and admin1:
            cent = centroids.us_state_centroid(admin1)
            if not cent:
                cent = centroids.country_centroid(country)
                key = (country, None)
                label = "USA"
            else:
                key = (country, admin1)
                label = f"USA — {admin1}"
        else:
            cent = centroids.country_centroid(country)
            key = (country, None)
            label = country

        if not cent:
            continue

        bucket = by_key.setdefault(key, {
            "label": label,
            "lon": cent[0],
            "lat": cent[1],
            "n_human": 0,
            "n_nonhuman": 0,
            "most_recent_iso": None,
        })
        if m.get("source_category") == "Human":
            bucket["n_human"] += 1
        else:
            bucket["n_nonhuman"] += 1

        # Track the most recent collection_date at this location.
        # ISO date strings compare correctly with string comparison.
        cd = m.get("collection_date")
        if cd:
            cd_str = str(cd)[:10]
            if not bucket["most_recent_iso"] or cd_str > bucket["most_recent_iso"]:
                bucket["most_recent_iso"] = cd_str

    # Compute days_ago from today for each location
    today_iso = today.isoformat()
    result = []
    for bucket in by_key.values():
        days_ago = None
        if bucket["most_recent_iso"]:
            try:
                most_recent = _date.fromisoformat(bucket["most_recent_iso"])
                days_ago = (today - most_recent).days
            except (ValueError, TypeError):
                days_ago = None
        bucket["days_ago"] = days_ago
        # Drop the raw iso string from the output — only days_ago is needed
        bucket.pop("most_recent_iso", None)
        result.append(bucket)

    return result
