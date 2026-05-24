"""Inline SVG histogram for cluster cards.

Renders a small stacked-bar histogram showing per-year human vs. nonhuman
isolate counts. Server-side, no JavaScript. Designed to be embedded
inline in HTML.

Color choice follows CDC-style epi curve convention:
  - Human (clinical) → deep red (rust)
  - Nonhuman (food/animal/environment/unknown) → teal

The two colors must be distinguishable to red-green colorblind viewers;
this red/teal pair passes that test.

The chart self-scales to its data: the y-axis tops out at max(year_total)
and the x-axis covers the full year range with one bar per year.
"""

from __future__ import annotations

from typing import Iterable


# Visual constants — keep small, refer to CSS variables in main.css where
# possible so the histogram inherits site theming.
W = 280                  # total svg width
H = 90                   # total svg height
PAD_LEFT = 28            # room for y-axis labels
PAD_RIGHT = 6
PAD_TOP = 10
PAD_BOTTOM = 22          # room for x-axis labels

COLOR_HUMAN = "#b91c1c"      # rust (matches --rust in main.css)
COLOR_NONHUMAN = "#0f766e"   # teal (matches --teal in main.css)
COLOR_AXIS = "#444"
COLOR_GRID = "#d8d2c4"


def render_histogram_svg(
    histogram: list[dict],
    max_year_count: int | None = None,
) -> str:
    """Render a two-color stacked bar histogram as an SVG string.

    `histogram` is a list of dicts: [{year: int, n_human: int, n_nonhuman: int}, ...]
    `max_year_count` is the max single-year total (for scaling); if None,
    computed from the data.
    """
    if not histogram:
        return _empty_svg()

    if max_year_count is None or max_year_count <= 0:
        max_year_count = max(
            (h["n_human"] + h["n_nonhuman"] for h in histogram),
            default=1,
        )
    if max_year_count <= 0:
        return _empty_svg()

    years = [h["year"] for h in histogram]
    year_min = min(years)
    year_max = max(years)
    # Always show a contiguous year range; fill missing years with zeros
    by_year = {h["year"]: h for h in histogram}
    rows = []
    for yr in range(year_min, year_max + 1):
        h = by_year.get(yr, {"year": yr, "n_human": 0, "n_nonhuman": 0})
        rows.append(h)

    n_bars = len(rows)
    chart_w = W - PAD_LEFT - PAD_RIGHT
    chart_h = H - PAD_TOP - PAD_BOTTOM
    bar_total_w = chart_w / n_bars
    bar_w = max(1.5, bar_total_w * 0.78)  # leave a small gap between bars
    bar_gap = bar_total_w - bar_w

    bars_svg: list[str] = []
    label_svg: list[str] = []

    for i, h in enumerate(rows):
        total = h["n_human"] + h["n_nonhuman"]
        x = PAD_LEFT + bar_gap / 2 + i * bar_total_w

        if total > 0:
            # Heights are proportional to count
            h_human = chart_h * (h["n_human"] / max_year_count)
            h_nonhuman = chart_h * (h["n_nonhuman"] / max_year_count)

            # Bars stack from bottom: nonhuman first (bottom), human on top
            y_nonhuman = PAD_TOP + chart_h - h_nonhuman
            y_human = y_nonhuman - h_human

            if h["n_nonhuman"] > 0:
                bars_svg.append(
                    f'<rect x="{x:.1f}" y="{y_nonhuman:.1f}" '
                    f'width="{bar_w:.1f}" height="{h_nonhuman:.1f}" '
                    f'fill="{COLOR_NONHUMAN}" />'
                )
            if h["n_human"] > 0:
                bars_svg.append(
                    f'<rect x="{x:.1f}" y="{y_human:.1f}" '
                    f'width="{bar_w:.1f}" height="{h_human:.1f}" '
                    f'fill="{COLOR_HUMAN}" />'
                )

        # X-axis year labels — show every year if ≤6 bars, every other if ≤12, then sparse.
        # Always show first and last.
        if (
            n_bars <= 6
            or (n_bars <= 12 and i % 2 == 0)
            or (n_bars <= 24 and i % 4 == 0)
            or (i % max(1, n_bars // 6) == 0)
            or i == n_bars - 1
        ):
            cx = x + bar_w / 2
            label_svg.append(
                f'<text x="{cx:.1f}" y="{H - 6}" font-size="9" '
                f'fill="{COLOR_AXIS}" text-anchor="middle">{h["year"]}</text>'
            )

    # Y-axis: a single tick label at the max
    y_axis_label = (
        f'<text x="{PAD_LEFT - 4}" y="{PAD_TOP + 8}" font-size="9" '
        f'fill="{COLOR_AXIS}" text-anchor="end">{max_year_count}</text>'
        f'<text x="{PAD_LEFT - 4}" y="{H - PAD_BOTTOM}" font-size="9" '
        f'fill="{COLOR_AXIS}" text-anchor="end">0</text>'
    )

    # Axes (lines)
    axes = (
        f'<line x1="{PAD_LEFT}" y1="{PAD_TOP}" x2="{PAD_LEFT}" y2="{H - PAD_BOTTOM}" '
        f'stroke="{COLOR_GRID}" stroke-width="1" />'
        f'<line x1="{PAD_LEFT}" y1="{H - PAD_BOTTOM}" x2="{W - PAD_RIGHT}" y2="{H - PAD_BOTTOM}" '
        f'stroke="{COLOR_AXIS}" stroke-width="1" />'
    )

    svg = (
        f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" '
        f'class="cluster-histogram" role="img" '
        f'aria-label="Yearly isolate counts: red=human, teal=nonhuman">'
        f'{axes}'
        f'{y_axis_label}'
        f'{"".join(bars_svg)}'
        f'{"".join(label_svg)}'
        f'</svg>'
    )
    return svg


def _empty_svg() -> str:
    """SVG placeholder when there's no histogram data."""
    return (
        f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" '
        f'class="cluster-histogram cluster-histogram-empty">'
        f'<text x="{W/2}" y="{H/2}" font-size="11" fill="{COLOR_AXIS}" '
        f'text-anchor="middle" dominant-baseline="middle" font-style="italic">'
        f'No dated isolates'
        f'</text></svg>'
    )
