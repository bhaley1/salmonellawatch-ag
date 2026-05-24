"""Listeria monocytogenes MLST → Clonal Complex lookup, and CC notes.

DRAFT — review before deploying to a public-facing site.

Clonal complexes (CCs) group together STs that descend from a common
ancestor (defined via goeBURST / eBURST clustering at the MLST level).
The PubMLST database publishes the canonical ST → CC assignments for
Listeria monocytogenes. This module reproduces the assignments for the
STs most commonly seen in outbreaks; STs not in this table return CC=None.

The "outbreak notes" below are intended as SUMMARY LABELS only, not as
clinical diagnostic claims. They reflect well-documented patterns from
the public-health literature but are deliberately phrased to avoid
deterministic statements.
"""

from __future__ import annotations


# ST (integer) → CC label. Sourced from PubMLST L. monocytogenes scheme.
# This is a curated subset focusing on STs that appear in outbreak literature;
# many other STs exist and would return CC=None here.
ST_TO_CC: dict[int, str] = {
    # CC1 — historically associated with severe disease
    1: "CC1", 3: "CC1",
    # CC2 — large multistate outbreaks
    2: "CC2", 145: "CC2",
    # CC4 — strong maternal-fetal association
    4: "CC4", 217: "CC4",
    # CC5
    5: "CC5",
    # CC6 — recently emerged, deli-meat outbreaks
    6: "CC6", 376: "CC6",
    # CC7
    7: "CC7",
    # CC8 — environmental persistence
    8: "CC8", 120: "CC8",
    # CC9 — abundant in food, lower virulence
    9: "CC9", 122: "CC9",
    # CC11
    11: "CC11",
    # CC14
    14: "CC14",
    # CC18
    18: "CC18",
    # CC19
    19: "CC19",
    # CC26
    26: "CC26",
    # CC31
    31: "CC31",
    # CC37
    37: "CC37",
    # CC54
    54: "CC54",
    # CC59
    59: "CC59",
    # CC87
    87: "CC87",
    # CC88
    88: "CC88",
    # CC101
    101: "CC101",
    # CC121 — persistent in food processing, recent clinical relevance
    121: "CC121",
    # CC155
    155: "CC155",
    # CC199
    199: "CC199",
    # CC204
    204: "CC204",
    # CC224
    224: "CC224",
    # CC315
    315: "CC315",
    # CC321
    321: "CC321",
    # CC382
    382: "CC382",
    # CC388
    388: "CC388",
    # CC403
    403: "CC403",
    # CC415
    415: "CC415",
    # CC499
    499: "CC499",
}


# CC → curated notes for the dashboard. Keep labels SHORT (will appear as
# inline chips) and the longer note as a tooltip. Phrase as observed
# patterns, not deterministic claims.
CC_NOTES: dict[str, dict[str, str]] = {
    "CC1": {
        "label": "hypervirulent lineage",
        "tooltip": "CC1 isolates historically over-represented in invasive listeriosis cases and associated with dairy-related outbreaks.",
    },
    "CC2": {
        "label": "outbreak-associated lineage",
        "tooltip": "CC2 has appeared in multiple large multistate outbreaks; clinically common.",
    },
    "CC4": {
        "label": "maternal-fetal lineage",
        "tooltip": "CC4 (the 'Listeria epidemica' lineage) shows strong association with maternal-fetal disease and severe outcomes.",
    },
    "CC6": {
        "label": "emerging outbreak lineage",
        "tooltip": "CC6 has been associated with recent deli-meat outbreaks; classified as hypervirulent in recent literature.",
    },
    "CC9": {
        "label": "food-common, lower-virulence",
        "tooltip": "CC9 is frequently isolated from food and food-processing environments but appears less often than expected in clinical cases.",
    },
    "CC121": {
        "label": "processing-environment persistent",
        "tooltip": "CC121 is well-documented as persistent in food-processing facilities over years to decades. Recent literature suggests increasing clinical relevance.",
    },
}


def st_to_cc(st_value: str | None) -> str | None:
    """Map an ST string (e.g. 'ST6') to a CC string (e.g. 'CC6') or None."""
    if not st_value:
        return None
    s = st_value.strip().upper()
    if s.startswith("ST"):
        s = s[2:]
    if not s.isdigit():
        return None
    return ST_TO_CC.get(int(s))


def cc_notes(cc: str | None) -> dict[str, str] | None:
    """Return {label, tooltip} for a CC, or None if no curated notes."""
    if not cc:
        return None
    return CC_NOTES.get(cc)
