"""Tests for the SVG map renderer."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.lookups import centroids
from src.render import map as cluster_map


# ---------------------------------------------------------------------------
# Centroid lookup
# ---------------------------------------------------------------------------

def test_country_centroid_known():
    c = centroids.country_centroid("United States")
    assert c is not None
    assert -180 <= c[0] <= 180
    assert -90 <= c[1] <= 90


def test_country_centroid_unknown():
    assert centroids.country_centroid("Ruritania") is None
    assert centroids.country_centroid(None) is None
    assert centroids.country_centroid("") is None


def test_us_state_centroid():
    assert centroids.us_state_centroid("Maryland") is not None
    assert centroids.us_state_centroid("Texas") is not None
    assert centroids.us_state_centroid("District of Columbia") is not None
    # Wrong: state shouldn't match
    assert centroids.us_state_centroid("Ontario") is None
    assert centroids.us_state_centroid(None) is None


# ---------------------------------------------------------------------------
# Projection
# ---------------------------------------------------------------------------

def test_projection_basic():
    """Equator/prime meridian should land in the middle of the canvas."""
    x, y = cluster_map._project(0, 0)
    # x ≈ middle of canvas
    assert abs(x - cluster_map.W / 2) < 2
    # y ≈ slightly south of middle (because LAT_MAX is bigger than -LAT_MIN)
    expected_y = (cluster_map.LAT_MAX - 0) / (cluster_map.LAT_MAX - cluster_map.LAT_MIN) * (cluster_map.H - 2 * cluster_map.PAD_Y) + cluster_map.PAD_Y
    assert abs(y - expected_y) < 1


def test_projection_extremes():
    """Corner projections should land in the corners."""
    x_west, _ = cluster_map._project(cluster_map.LON_MIN, 0)
    x_east, _ = cluster_map._project(cluster_map.LON_MAX, 0)
    assert x_west < x_east
    assert abs(x_west - cluster_map.PAD_X) < 0.1
    assert abs(x_east - (cluster_map.W - cluster_map.PAD_X)) < 0.1


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def test_aggregate_us_by_state():
    members = [
        {"geo_country": "United States", "geo_admin1": "Maryland",
         "source_category": "Human"},
        {"geo_country": "United States", "geo_admin1": "Texas",
         "source_category": "Food"},
    ]
    locs = cluster_map.aggregate_locations(members)
    labels = {l["label"] for l in locs}
    assert "USA — Maryland" in labels
    assert "USA — Texas" in labels


def test_aggregate_us_unknown_state_falls_back():
    """A US isolate with an unrecognized admin1 should still place at country level."""
    members = [
        {"geo_country": "United States", "geo_admin1": "Foovania",
         "source_category": "Human"},
    ]
    locs = cluster_map.aggregate_locations(members)
    assert len(locs) == 1
    assert locs[0]["label"] == "USA"


def test_aggregate_non_us_by_country():
    members = [
        {"geo_country": "Germany", "geo_admin1": "Bavaria",
         "source_category": "Human"},
        {"geo_country": "Germany", "geo_admin1": None,
         "source_category": "Food"},
    ]
    locs = cluster_map.aggregate_locations(members)
    # Both rolled up to one country-level dot
    assert len(locs) == 1
    assert locs[0]["label"] == "Germany"
    assert locs[0]["n_human"] == 1
    assert locs[0]["n_nonhuman"] == 1


def test_aggregate_split_counts():
    """Human/nonhuman counts should split correctly at each location."""
    members = [
        {"geo_country": "United States", "geo_admin1": "Maryland", "source_category": "Human"},
        {"geo_country": "United States", "geo_admin1": "Maryland", "source_category": "Human"},
        {"geo_country": "United States", "geo_admin1": "Maryland", "source_category": "Food"},
        {"geo_country": "United States", "geo_admin1": "Texas", "source_category": "Environment"},
    ]
    locs = cluster_map.aggregate_locations(members)
    by_label = {l["label"]: l for l in locs}
    assert by_label["USA — Maryland"]["n_human"] == 2
    assert by_label["USA — Maryland"]["n_nonhuman"] == 1
    assert by_label["USA — Texas"]["n_human"] == 0
    assert by_label["USA — Texas"]["n_nonhuman"] == 1


def test_aggregate_drops_unmappable():
    """Members without a geocodable country are silently dropped."""
    members = [
        {"geo_country": "Ruritania", "geo_admin1": None, "source_category": "Human"},
        {"geo_country": "United States", "geo_admin1": "Maryland", "source_category": "Human"},
    ]
    locs = cluster_map.aggregate_locations(members)
    assert len(locs) == 1
    assert locs[0]["label"] == "USA — Maryland"


# ---------------------------------------------------------------------------
# SVG rendering
# ---------------------------------------------------------------------------

def test_render_empty():
    """Empty locations should render the placeholder SVG."""
    svg = cluster_map.render_cluster_map_svg([])
    assert svg.startswith("<svg")
    assert "No mapped isolates" in svg


def test_opacity_recency():
    """Dots fade based on most-recent isolate age."""
    # Recent (within 1 year): full opacity
    assert cluster_map._opacity_for(0) == 1.0
    assert cluster_map._opacity_for(365) == 1.0
    # Old (over 5 years): minimum opacity
    assert cluster_map._opacity_for(1825) == 0.35
    assert cluster_map._opacity_for(10000) == 0.35
    # Middle: interpolated
    mid = cluster_map._opacity_for(1095)  # ~3 years
    assert 0.35 < mid < 1.0
    # Unknown: medium fade
    assert cluster_map._opacity_for(None) == 0.55


def test_aggregate_tracks_days_ago():
    """Aggregator records days since most-recent collection per location."""
    from datetime import date
    today = date(2026, 5, 22)
    members = [
        {"geo_country": "United States", "geo_admin1": "Maryland",
         "source_category": "Human", "collection_date": "2026-05-15"},
        {"geo_country": "United States", "geo_admin1": "Maryland",
         "source_category": "Food", "collection_date": "2020-03-01"},
        {"geo_country": "Germany", "geo_admin1": None,
         "source_category": "Human", "collection_date": "2010-01-01"},
    ]
    locs = cluster_map.aggregate_locations(members, today=today)
    by_label = {l["label"]: l for l in locs}
    # Maryland: most recent is 2026-05-15 → days_ago = 7
    assert by_label["USA — Maryland"]["days_ago"] == 7
    # Germany: most recent is 2010-01-01 → many years ago
    assert by_label["Germany"]["days_ago"] > 5000


def test_aggregate_no_date_yields_none():
    """A member without collection_date doesn't crash; days_ago becomes None."""
    members = [
        {"geo_country": "United States", "geo_admin1": "Maryland",
         "source_category": "Human", "collection_date": None},
    ]
    locs = cluster_map.aggregate_locations(members)
    assert locs[0]["days_ago"] is None


def test_render_applies_opacity_to_old_dots():
    """An old location should render with reduced opacity."""
    locations = [
        {"label": "USA — Maryland", "lon": -76.8, "lat": 39.0,
         "n_human": 3, "n_nonhuman": 0, "days_ago": 10},      # recent → 1.0
        {"label": "Germany", "lon": 10.4, "lat": 51.2,
         "n_human": 5, "n_nonhuman": 5, "days_ago": 3000},     # old → faded
    ]
    svg = cluster_map.render_cluster_map_svg(locations)
    # The German dot should carry an opacity attribute < 1
    assert 'opacity="' in svg


def test_render_with_data():
    locations = [
        {"label": "USA — Maryland", "lon": -76.8, "lat": 39.0,
         "n_human": 3, "n_nonhuman": 0},
        {"label": "USA — Texas", "lon": -99.3, "lat": 31.0,
         "n_human": 0, "n_nonhuman": 5},
        {"label": "Germany", "lon": 10.4, "lat": 51.2,
         "n_human": 2, "n_nonhuman": 4},
    ]
    svg = cluster_map.render_cluster_map_svg(locations)
    assert "<svg" in svg
    # Human-only dot uses HUMAN color
    assert cluster_map.COLOR_HUMAN in svg
    # Nonhuman-only and mixed both use NONHUMAN color
    assert cluster_map.COLOR_NONHUMAN in svg
    # Tooltips should include location labels
    assert "USA — Maryland: 3 human" in svg
    assert "USA — Texas: 0 human" in svg
    assert "Germany: 2 human" in svg


def test_render_single_location():
    """One location should render cleanly with one dot."""
    locations = [
        {"label": "Iceland", "lon": -18.9, "lat": 64.9,
         "n_human": 1, "n_nonhuman": 0},
    ]
    svg = cluster_map.render_cluster_map_svg(locations)
    assert "Iceland: 1 human" in svg


def test_render_mixed_dot_has_both_colors():
    """A location with both human and nonhuman should produce a split dot."""
    locations = [
        {"label": "USA — Maryland", "lon": -76.8, "lat": 39.0,
         "n_human": 3, "n_nonhuman": 2},
    ]
    svg = cluster_map.render_cluster_map_svg(locations)
    # Mixed dot uses both colors in path fills
    assert f'fill="{cluster_map.COLOR_HUMAN}"' in svg
    assert f'fill="{cluster_map.COLOR_NONHUMAN}"' in svg
    assert "map-dot-mixed" in svg


def test_radius_scaling():
    """Bigger counts produce bigger dots."""
    r_small = cluster_map._radius_for(1)
    r_med = cluster_map._radius_for(10)
    r_big = cluster_map._radius_for(100)
    assert r_small < r_med < r_big
    # Clamped at MAX_R
    assert cluster_map._radius_for(10000) == cluster_map.MAX_R


# ---------------------------------------------------------------------------
# Test runner
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import inspect
    mod = sys.modules[__name__]
    tests = [
        (name, fn) for name, fn in inspect.getmembers(mod, inspect.isfunction)
        if name.startswith("test_")
    ]
    passed, failed = 0, 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
            passed += 1
        except Exception as e:
            import traceback
            print(f"  FAIL  {name}: {type(e).__name__}: {e}")
            traceback.print_exc()
            failed += 1
    print()
    print(f"{passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
