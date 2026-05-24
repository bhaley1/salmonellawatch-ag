"""Parse NCBI's TSV files into typed records ready for upsert.

This module owns:
  - Reading metadata.tsv into Isolate records
  - Reading cluster_list.tsv into PDT→PDS mappings
  - Reading amr.metadata.tsv into per-isolate gene calls
  - Source category inference (Human / Food / Animal / Environment / Unknown)
  - Country and admin-1 parsing from geo_loc_name
  - Date parsing with preservation of original strings
"""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
from typing import Iterator

from .. import config

log = logging.getLogger(__name__)


@dataclass
class IsolateRow:
    """A parsed isolate, ready for SQLite upsert. Mirrors the schema."""
    pdt_acc: str
    pathogen: str
    pds_acc: str | None
    epi_type: str | None
    source_category: str
    host: str | None
    isolation_source: str | None
    geo_loc_name: str | None
    geo_country: str | None
    geo_admin1: str | None
    collection_date: date | None
    collection_date_raw: str | None
    target_creation_date: date | None
    scientific_name: str | None
    serovar: str | None
    biosample_acc: str | None
    asm_acc: str | None
    sra_acc: str | None
    asm_level: str | None
    asm_stats_contig_n50: int | None
    food_origin: str | None
    ifsac_category: str | None
    host_disease: str | None
    bioproject_acc: str | None
    pdg_release: str


@dataclass
class AmrRow:
    pdt_acc: str
    gene_symbol: str
    element_subtype: str | None


# ---------------------------------------------------------------------------
# Column lookup — NCBI's schema has drifted over time; keep candidates
# ---------------------------------------------------------------------------
COL_CANDIDATES: dict[str, tuple[str, ...]] = {
    "pdt_acc":              ("target_acc", "PDT_acc", "pdt_acc"),
    "epi_type":             ("epi_type",),
    "host":                 ("host",),
    "isolation_source":     ("isolation_source", "source"),
    "geo_loc_name":         ("geo_loc_name", "country", "location"),
    "collection_date":      ("collection_date",),
    "target_creation_date": ("target_creation_date",),
    "scientific_name":      ("scientific_name", "organism"),
    "serovar":              ("serovar", "serotype"),
    "biosample_acc":        ("biosample_acc", "BioSample"),
    "asm_acc":              ("asm_acc",),
    "sra_acc":              ("Run", "sra_acc"),
    "asm_level":            ("asm_level",),
    "asm_stats_contig_n50": ("asm_stats_contig_n50",),
    "food_origin":          ("food_origin",),
    "ifsac_category":       ("IFSAC_category", "ifsac_category"),
    "host_disease":         ("host_disease",),
    "bioproject_acc":       ("bioproject_acc", "BioProject"),
}


def _first_present(row: dict, names: tuple[str, ...]) -> str | None:
    for n in names:
        v = row.get(n)
        if v is not None and v != "" and str(v).upper() != "NULL":
            return v.strip() if isinstance(v, str) else v
    return None


def _to_date(s: str | None) -> date | None:
    if not s or str(s).strip() == "" or str(s).upper() == "NULL":
        return None
    s = str(s).strip()
    if "T" in s:
        s = s.split("T", 1)[0]
    if " " in s:
        s = s.split(" ", 1)[0]
    for fmt in ("%Y-%m-%d", "%Y-%m", "%Y/%m/%d", "%Y"):
        try:
            return datetime.strptime(s, fmt).date()
        except ValueError:
            continue
    return None


def _to_int(v) -> int | None:
    if v is None or v == "" or str(v).upper() == "NULL":
        return None
    try:
        return int(v)
    except (ValueError, TypeError):
        return None


# ---------------------------------------------------------------------------
# Geography parsing
# ---------------------------------------------------------------------------
# NCBI's geo_loc_name is typically "Country: Admin1" or "Country: Admin1, Admin2"
# or just "Country". We parse out country and (where available) admin1.

def _parse_geo(geo: str | None) -> tuple[str | None, str | None]:
    if not geo:
        return (None, None)
    s = geo.strip()
    if ":" in s:
        country, rest = s.split(":", 1)
        country = country.strip() or None
        admin1 = rest.split(",")[0].strip() if rest else None
        admin1 = admin1 or None
        return (country, admin1)
    return (s, None)


# ---------------------------------------------------------------------------
# Source category inference
# ---------------------------------------------------------------------------
# We classify into 5 buckets: Human, Food, Animal, Environment, Unknown.
# The priority is:
#   1. epi_type == 'clinical' → Human
#   2. host matches a human token → Human
#   3. isolation_source matches a human-source token → Human
#   4. isolation_source matches food tokens → Food
#   5. host or isolation_source matches animal tokens → Animal
#   6. isolation_source matches environment tokens → Environment
#   7. Otherwise → Unknown
#
# These heuristics are conservative; "Unknown" is the right answer when
# NCBI's metadata is genuinely ambiguous.

def _matches_any(text: str | None, tokens: set[str]) -> bool:
    if not text:
        return False
    t = text.lower()
    return any(tok in t for tok in tokens)


def infer_source_category(
    epi_type: str | None,
    host: str | None,
    isolation_source: str | None,
) -> str:
    """Classify an isolate's source for surveillance summarization.

    Priority order is deliberate:
      1. Human — clinical samples and human hosts always win
      2. Environment — facility/water/soil swabs are a distinct surveillance
         signal from food products, even if the words "food" or "milk" appear
         in the context (e.g., "food processing facility drain")
      3. Food — explicit food product terms
      4. Animal — animal hosts without food/environment context
      5. Unknown — when nothing matches confidently

    A subtle wrinkle: words like "feces", "blood", "urine" are in
    HUMAN_SOURCE_TOKENS because they're clinical samples, but they also
    apply to animals. So if the host is clearly an animal, don't classify
    the isolate as Human even if the source string contains those terms.
    """
    if (epi_type or "").lower() == "clinical":
        return "Human"
    if _matches_any(host, config.HUMAN_HOST_TOKENS):
        return "Human"
    host_is_animal = _matches_any(host, config.ANIMAL_TOKENS)
    if not host_is_animal and _matches_any(isolation_source, config.HUMAN_SOURCE_TOKENS):
        return "Human"
    # Environment before Food: a "food processing facility drain" is
    # facility surveillance, not a food sample. The Environment tokens
    # ("processing", "facility", "drain", "swab", "wastewater") are
    # specific enough that they should dominate over generic food words.
    if _matches_any(isolation_source, config.ENVIRONMENT_TOKENS):
        return "Environment"
    if _matches_any(isolation_source, config.FOOD_TOKENS):
        return "Food"
    if host_is_animal or _matches_any(isolation_source, config.ANIMAL_TOKENS):
        return "Animal"
    return "Unknown"


def infer_epi_type(host: str | None, isolation_source: str | None) -> str | None:
    """Best-effort epi_type if the column wasn't populated."""
    if _matches_any(host, config.HUMAN_HOST_TOKENS):
        return "clinical"
    if _matches_any(isolation_source, config.HUMAN_SOURCE_TOKENS):
        return "clinical"
    if _matches_any(isolation_source, config.FOOD_TOKENS):
        return "environmental/other"
    if _matches_any(host, config.ANIMAL_TOKENS) or _matches_any(isolation_source, config.ANIMAL_TOKENS):
        return "environmental/other"
    if _matches_any(isolation_source, config.ENVIRONMENT_TOKENS):
        return "environmental/other"
    return None


# ---------------------------------------------------------------------------
# cluster_list parser — PDT_acc → PDS_acc
# ---------------------------------------------------------------------------
def load_pds_map(path: Path | None) -> dict[str, str]:
    if not path or not path.exists():
        return {}
    out: dict[str, str] = {}
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fns = reader.fieldnames or []
        pds_col = next(
            (c for c in fns if c.lower() in ("pds_acc", "snp_cluster", "cluster")),
            None,
        )
        pdt_col = next(
            (c for c in fns if c.lower() in ("target_acc", "pdt_acc")),
            None,
        )
        if not pds_col or not pdt_col:
            log.warning("cluster_list missing PDS/PDT columns; got %s", fns)
            return {}
        for row in reader:
            pdt = (row.get(pdt_col) or "").strip()
            pds = (row.get(pds_col) or "").strip()
            if pdt and pds and pds.upper() != "NULL":
                out[pdt] = pds
    log.info("Loaded %d PDT→PDS mappings from %s", len(out), path.name)
    return out


# ---------------------------------------------------------------------------
# AMR parser
# ---------------------------------------------------------------------------
def iter_amr_rows(path: Path | None) -> Iterator[AmrRow]:
    if not path or not path.exists():
        return
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fns = reader.fieldnames or []
        pdt_col = next((c for c in fns if c.lower() in ("target_acc", "pdt_acc")), None)
        sym_col = next(
            (c for c in fns if c.lower() in (
                "element_symbol", "gene_symbol", "symbol", "amr_gene", "gene",
            )),
            None,
        )
        subtype_col = next(
            (c for c in fns if c.lower() in ("element_subtype", "subtype", "type")),
            None,
        )
        if not pdt_col or not sym_col:
            log.warning("AMR file missing key columns; got %s", fns)
            return
        for row in reader:
            pdt = (row.get(pdt_col) or "").strip()
            sym = (row.get(sym_col) or "").strip()
            if not pdt or not sym:
                continue
            subtype = (row.get(subtype_col) or "").strip() if subtype_col else None
            yield AmrRow(pdt_acc=pdt, gene_symbol=sym, element_subtype=subtype or None)


# ---------------------------------------------------------------------------
# metadata.tsv parser
# ---------------------------------------------------------------------------
def iter_isolates(
    metadata_path: Path,
    pathogen: str,
    pds_map: dict[str, str],
    pdg_release: str,
) -> Iterator[IsolateRow]:
    """Yield IsolateRow objects parsed from metadata.tsv."""
    with open(metadata_path, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        log.info(
            "[%s] TSV columns: %s",
            pathogen,
            ", ".join((reader.fieldnames or [])[:60]),
        )
        for row in reader:
            pdt = _first_present(row, COL_CANDIDATES["pdt_acc"])
            if not pdt:
                continue

            host = _first_present(row, COL_CANDIDATES["host"])
            iso_src = _first_present(row, COL_CANDIDATES["isolation_source"])
            epi_raw = _first_present(row, COL_CANDIDATES["epi_type"])
            epi = epi_raw if epi_raw else infer_epi_type(host, iso_src)
            source_category = infer_source_category(epi, host, iso_src)

            geo_raw = _first_present(row, COL_CANDIDATES["geo_loc_name"])
            country, admin1 = _parse_geo(geo_raw)

            coll_raw = _first_present(row, COL_CANDIDATES["collection_date"])
            tgt_creation_raw = _first_present(row, COL_CANDIDATES["target_creation_date"])

            yield IsolateRow(
                pdt_acc=pdt,
                pathogen=pathogen,
                pds_acc=pds_map.get(pdt),
                epi_type=epi,
                source_category=source_category,
                host=host,
                isolation_source=iso_src,
                geo_loc_name=geo_raw,
                geo_country=country,
                geo_admin1=admin1,
                collection_date=_to_date(coll_raw),
                collection_date_raw=coll_raw,
                target_creation_date=_to_date(tgt_creation_raw),
                scientific_name=_first_present(row, COL_CANDIDATES["scientific_name"]),
                serovar=_first_present(row, COL_CANDIDATES["serovar"]),
                biosample_acc=_first_present(row, COL_CANDIDATES["biosample_acc"]),
                asm_acc=_first_present(row, COL_CANDIDATES["asm_acc"]),
                sra_acc=_first_present(row, COL_CANDIDATES["sra_acc"]),
                asm_level=_first_present(row, COL_CANDIDATES["asm_level"]),
                asm_stats_contig_n50=_to_int(
                    _first_present(row, COL_CANDIDATES["asm_stats_contig_n50"])
                ),
                food_origin=_first_present(row, COL_CANDIDATES["food_origin"]),
                ifsac_category=_first_present(row, COL_CANDIDATES["ifsac_category"]),
                host_disease=_first_present(row, COL_CANDIDATES["host_disease"]),
                bioproject_acc=_first_present(row, COL_CANDIDATES["bioproject_acc"]),
                pdg_release=pdg_release,
            )
