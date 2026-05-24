"""Clinically important AMR gene flags for Listeria monocytogenes.

DRAFT — review before deploying to a public-facing site.

Listeria is intrinsically resistant to several drug classes (cephalosporins,
fosfomycin, fluoroquinolones to varying degrees). What clinicians care about
is ACQUIRED resistance to first-line treatments:

  First-line invasive listeriosis:    ampicillin or penicillin + gentamicin
  Penicillin-allergic alternative:    trimethoprim-sulfamethoxazole

Any acquired resistance to those classes is a real clinical concern.
Listeria acquired resistance is uncommon overall, so flagging is rare and
therefore meaningful when it appears.

This module flags clusters that carry a "first-line resistance gene."
That phrasing is deliberately conservative: the flag means a resistance
gene is present in NCBI's AMRFinderPlus output, NOT that the strain is
phenotypically resistant in a clinical microbiology lab.

Gene symbols match AMRFinderPlus naming conventions as of mid-2026.
"""

from __future__ import annotations

# Class → set of gene symbols (or symbol prefixes ending in *)
FIRST_LINE_RESISTANCE_GENES: dict[str, set[str]] = {
    "Tetracyclines": {
        "tet(M)", "tet(O)", "tet(S)", "tet(T)", "tet(W)", "tet(L)", "tet(K)",
    },
    "Macrolides / lincosamides": {
        "erm(A)", "erm(B)", "erm(C)", "erm(G)", "mef(A)", "mef(E)",
        "lnu(A)", "lnu(B)", "lnu(G)",
    },
    "Aminoglycosides": {
        "aph(3')-III", "aph(3')-Ia", "aph(2'')-Ia",
        "aac(6')-aph(2'')", "aac(6')-Ie-aph(2'')-Ia",
        "aadA", "aadE", "ant(6)-Ia",
    },
    "Trimethoprim / sulfonamides": {
        "dfrG", "dfrD", "dfrK", "dfrA",
        "sul1", "sul2", "sul3",
    },
    "Chloramphenicol / phenicols": {
        "cat", "catA", "catB", "fexA", "fexB",
    },
    "Penicillins (rare in Listeria)": {
        # Beta-lactamases of concern if seen in Listeria
        "blaZ", "mecA", "mecC",
    },
}


def amr_clinical_classes(gene_symbols: set[str]) -> list[str]:
    """Return the set of clinical drug classes for which the input gene
    set contains a first-line resistance marker."""
    classes_present: list[str] = []
    # Normalize input: handle "tet(M)" variants like "tetM" stripped of parens
    normalized_in = {g.strip() for g in gene_symbols}
    normalized_lower = {g.lower() for g in normalized_in}
    for cls, genes in FIRST_LINE_RESISTANCE_GENES.items():
        for gene in genes:
            if (
                gene in normalized_in
                or gene.lower() in normalized_lower
            ):
                classes_present.append(cls)
                break
    return classes_present
