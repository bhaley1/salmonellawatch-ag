"""Salmonella enterica MLST → lineage/clade lookup and outbreak notes.

Salmonella does not use clonal complexes (CCs) in the same way as Listeria.
Instead, the epidemiologically meaningful groupings are serotype + ST.
This module maps common outbreak-associated STs to short descriptive labels
for display on the dashboard.

ST assignments sourced from PubMLST senterica scheme and outbreak literature.
STs not in this table return CC=None (displayed as ST only, no label).
"""

from __future__ import annotations


# ST (integer) → short lineage label
ST_TO_CC: dict[int, str] = {
    # Enteritidis — most common human salmonellosis globally
    11: "ST11",   # S. Enteritidis dominant ST
    # Typhimurium — second most common, broad host range
    19: "ST19",   # S. Typhimurium dominant ST
    34: "ST34",   # S. Typhimurium variant, monophasic, high AMR
    # Newport — common in cattle, produce outbreaks
    5: "ST5",     # S. Newport
    # Heidelberg — poultry-associated, frequent MDR
    15: "ST15",   # S. Heidelberg
    # Infantis — broiler-associated, emerging MDR globally
    32: "ST32",   # S. Infantis
    # Kentucky — poultry, high ciprofloxacin resistance
    198: "ST198",  # S. Kentucky ciprofloxacin-resistant lineage
    # Dublin — cattle, invasive disease, high mortality
    10: "ST10",   # S. Dublin
    # Javiana — produce/amphibian-associated, SE USA
    167: "ST167",
    # Bareilly — seafood-associated
    182: "ST182",
    # Thompson — seafood/produce
    26: "ST26",
    # Oranienburg — sugar/confectionery outbreaks
    44: "ST44",
    # Stanley — peanut butter outbreaks
    1751: "ST1751",
    # Reading — baby food outbreaks Europe
    1625: "ST1625",
}

# Outbreak notes for display tooltips
CC_NOTES: dict[str, str] = {
    "ST11":   "S. Enteritidis — dominant global lineage; poultry and eggs",
    "ST19":   "S. Typhimurium — broad host range; beef, pork, produce",
    "ST34":   "S. Typhimurium monophasic variant — high AMR, emerging globally",
    "ST5":    "S. Newport — cattle and produce; MDR strains common",
    "ST15":   "S. Heidelberg — poultry-associated; frequent MDR",
    "ST32":   "S. Infantis — broiler chicken; pESI-like megaplasmid MDR",
    "ST198":  "S. Kentucky — ciprofloxacin-resistant lineage; poultry",
    "ST10":   "S. Dublin — cattle; invasive disease; high case-fatality",
    "ST167":  "S. Javiana — produce and amphibians; SE United States",
    "ST182":  "S. Bareilly — seafood-associated",
    "ST26":   "S. Thompson — seafood and produce",
    "ST44":   "S. Oranienburg — confectionery and produce",
    "ST1751": "S. Stanley — peanut butter and tree nut outbreaks",
    "ST1625": "S. Reading — infant formula; European outbreaks",
}


def st_to_cc(st: int | None) -> str | None:
    """Return the lineage label for a given ST, or None."""
    if st is None:
        return None
    return ST_TO_CC.get(st)


def cc_note(cc: str | None) -> str | None:
    """Return outbreak/epidemiology note for a CC/lineage label."""
    if cc is None:
        return None
    return CC_NOTES.get(cc)
