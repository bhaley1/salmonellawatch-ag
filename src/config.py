"""Configuration for Pathogen Watch v2.

Single source of truth for paths, pathogen definitions, and tuning knobs.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = REPO_ROOT / "src"
TEMPLATES_DIR = REPO_ROOT / "templates"
DB_DIR = REPO_ROOT / "db"
DB_PATH = DB_DIR / "pathogen-watch.sqlite"
SCHEMA_PATH = SRC_DIR / "schema.sql"
CACHE_DIR = REPO_ROOT / "cache"
SITE_DIR = REPO_ROOT / "site"


# ---------------------------------------------------------------------------
# NCBI Pathogen Detection
# ---------------------------------------------------------------------------
NCBI_BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results"


@dataclass(frozen=True)
class PathogenConfig:
    """Per-pathogen configuration. `taxgroup` is the NCBI directory name."""
    display_name: str        # 'Listeria', shown in UI
    taxgroup: str            # 'Listeria', NCBI's directory name
    scientific_name: str     # 'Listeria monocytogenes'
    stx_filter: bool = False # only True for STEC
    snps_subdir: str | None = None  # e.g. 'latest_snps' for Salmonella


PATHOGENS: dict[str, PathogenConfig] = {
    "Salmonella": PathogenConfig(
        display_name="Salmonella",
        taxgroup="Salmonella",
        scientific_name="Salmonella enterica",
        snps_subdir="latest_snps",
    ),
}


# ---------------------------------------------------------------------------
# Surveillance window
# ---------------------------------------------------------------------------
# 60 days reflects Salmonella's epidemiology: a 6-72 hour incubation
# period with typical submission lag of 2-4 weeks from state public health
# labs to NCBI. Cases visible in the last 60 days of target_creation_date
# correspond to exposures roughly 3-8 weeks ago. CDC uses 60-day windows
# for Salmonella outbreak investigation.
RECENT_WINDOW_DAYS = 180

# How far back the per-cluster geographic-footprint map looks. Longer than
# the recent-cases window because the map shows the cluster's broader
# recent geography, not just the cases that drove it onto the dashboard.
MAP_WINDOW_DAYS = 365


# ---------------------------------------------------------------------------
# Cache freshness
# ---------------------------------------------------------------------------
CACHE_MAX_AGE_HOURS = 3


# ---------------------------------------------------------------------------
# Review mode banner
# ---------------------------------------------------------------------------
# When True, the dashboard shows a prominent "PRIVATE REVIEW — DO NOT
# SHARE" banner. Use this while sharing with colleagues for feedback
# before public release. Toggle to False when ready to go public.
#
# Override via environment variable LW_REVIEW_MODE=0 or =1.
REVIEW_MODE = os.environ.get("LW_REVIEW_MODE", "1") == "1"


# ---------------------------------------------------------------------------
# Source category inference
# ---------------------------------------------------------------------------
# Used to bucket isolation_source/host strings into broad categories
# for surveillance summaries. These are heuristics, deliberately
# conservative — unknown sources stay 'Unknown' rather than being guessed.

HUMAN_HOST_TOKENS = {
    "homo sapiens", "human", "hu", "patient", "h. sapiens",
}
HUMAN_SOURCE_TOKENS = {
    "stool", "feces", "faeces", "blood", "csf", "urine", "wound",
    "sputum", "respiratory", "throat", "clinical", "patient",
}
FOOD_TOKENS = {
    "beef", "poultry", "chicken", "pork", "swine pork", "turkey",
    "deli", "meat", "ground beef", "ground pork", "ground turkey",
    "dairy", "milk", "cheese", "yogurt", "butter", "cream",
    "egg", "eggs",
    "produce", "lettuce", "spinach", "salad", "sprout", "sprouts",
    "fruit", "vegetable", "vegetables",
    "fish", "seafood", "shrimp", "salmon", "tuna", "oyster",
    "raw meat", "ready-to-eat",
    "deli meat", "salami", "pepperoni", "hot dog", "frankfurter",
    "ice cream",
    "food",
}
ANIMAL_TOKENS = {
    "cattle", "bovine", "cow", "calf", "heifer", "steer",
    "swine", "pig", "piglet", "sow", "boar",
    "poultry", "chicken", "broiler", "layer", "hen", "rooster",
    "turkey",  # may also indicate food product; food tokens match first
    "sheep", "lamb", "ovine",
    "goat", "caprine",
    "horse", "equine",
    "dog", "canine",
    "cat", "feline",
    "deer", "elk", "wildlife", "wild bird",
}
ENVIRONMENT_TOKENS = {
    "water", "soil", "sediment", "environmental", "environment",
    "swab", "drain", "floor", "wall", "ceiling",
    "compost", "manure", "wastewater", "sewage",
    "surface", "processing", "facility", "plant", "factory",
    "irrigation", "stream", "river", "lake", "ocean",
}
