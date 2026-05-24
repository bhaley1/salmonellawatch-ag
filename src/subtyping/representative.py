"""Select a representative isolate for MLST typing.

For each cluster we need to pick one isolate to run MLST against. The
chosen isolate determines the cluster's ST call (which is then propagated
to all members). Selection criteria, in priority order:

  1. Has a non-null asm_acc (NCBI assembly accession we can download)
  2. Highest asm_level: Complete Genome > Chromosome > Scaffold > Contig
  3. Highest asm_stats_contig_n50 (larger contigs = better assembly)
  4. Tiebreak: lowest pdt_acc lexically (stable, deterministic)

Returns None if no member has a downloadable assembly — those clusters
get mlst_error='no_assembly' in cluster_typing.
"""

from __future__ import annotations

import logging
import sqlite3
from dataclasses import dataclass

log = logging.getLogger(__name__)


# Higher is better — used to sort representatives.
_ASM_LEVEL_RANK = {
    "Complete Genome": 4,
    "Chromosome":      3,
    "Scaffold":        2,
    "Contig":          1,
}


@dataclass
class Representative:
    pathogen: str
    pds_acc: str
    pdt_acc: str
    asm_acc: str
    asm_level: str | None
    n50: int | None


def _sort_key(row) -> tuple:
    """Sort representatives best-first."""
    asm_level = row["asm_level"] or ""
    level_rank = _ASM_LEVEL_RANK.get(asm_level, 0)
    n50 = row["asm_stats_contig_n50"] or 0
    pdt = row["pdt_acc"] or ""
    # Negative for desc sort on the first two fields, positive on the tiebreak
    return (-level_rank, -n50, pdt)


def pick_representative(
    conn: sqlite3.Connection,
    pathogen: str,
    pds_acc: str,
) -> Representative | None:
    """Pick the best assembly-bearing isolate for a cluster, or None."""
    rows = conn.execute("""
        SELECT pdt_acc, asm_acc, asm_level, asm_stats_contig_n50
        FROM isolates
        WHERE pathogen = ? AND pds_acc = ?
              AND asm_acc IS NOT NULL AND asm_acc != ''
    """, (pathogen, pds_acc)).fetchall()

    if not rows:
        return None

    rows = sorted(rows, key=_sort_key)
    best = rows[0]
    return Representative(
        pathogen=pathogen,
        pds_acc=pds_acc,
        pdt_acc=best["pdt_acc"],
        asm_acc=best["asm_acc"],
        asm_level=best["asm_level"],
        n50=best["asm_stats_contig_n50"],
    )


def list_clusters_needing_typing(
    conn: sqlite3.Connection,
    pathogen: str | None = None,
    limit: int | None = None,
) -> list[tuple[str, str]]:
    """List (pathogen, pds_acc) tuples that don't yet have MLST typing.

    Clusters with mlst_error set are NOT retyped (failure is sticky to
    avoid infinite retries on assemblies that just don't work).
    """
    where = "cs.n_total > 0 AND (ct.mlst_st IS NULL AND ct.mlst_error IS NULL)"
    params: list = []
    if pathogen:
        where += " AND cs.pathogen = ?"
        params.append(pathogen)

    query = f"""
        SELECT cs.pathogen, cs.pds_acc
        FROM cluster_summary cs
        LEFT JOIN cluster_typing ct
               ON cs.pathogen = ct.pathogen AND cs.pds_acc = ct.pds_acc
        WHERE {where}
        ORDER BY cs.n_human DESC, cs.n_total DESC
    """
    if limit:
        query += " LIMIT ?"
        params.append(limit)

    rows = conn.execute(query, params).fetchall()
    return [(r["pathogen"], r["pds_acc"]) for r in rows]
