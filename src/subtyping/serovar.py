"""Compute cluster-level consensus serovar.

The serovar/serotype field is per-isolate in NCBI's metadata, but a
cluster is by definition a group of genomically related isolates that
should share evolutionary history. So one representative serovar should
apply to all members within a SNP threshold.

We compute the *modal* (most common) non-null serovar across the cluster
and record both the count that agreed and the total that had any serovar
value. The renderer can display confidence as "1/2a (245 of 267 agree)".
"""

from __future__ import annotations

import logging
import sqlite3
from collections import Counter
from dataclasses import dataclass

log = logging.getLogger(__name__)


@dataclass
class SerovarConsensus:
    """Result of consensus serovar computation for one cluster."""
    consensus_serovar: str | None      # None if no cluster member had a serovar
    n_agreed: int                       # how many members had the modal value
    n_total_with_serovar: int           # how many members had any non-null serovar

    @property
    def agreement_fraction(self) -> float | None:
        if self.n_total_with_serovar == 0:
            return None
        return self.n_agreed / self.n_total_with_serovar


def compute_consensus_serovar(serovars: list[str | None]) -> SerovarConsensus:
    """Compute consensus from a list of per-isolate serovar values."""
    populated = [s for s in serovars if s and s.strip()]
    if not populated:
        return SerovarConsensus(None, 0, 0)
    counter = Counter(populated)
    modal, n = counter.most_common(1)[0]
    return SerovarConsensus(
        consensus_serovar=modal,
        n_agreed=n,
        n_total_with_serovar=len(populated),
    )


def compute_all_cluster_serovars(
    conn: sqlite3.Connection,
) -> list[tuple[str, str, SerovarConsensus]]:
    """For every (pathogen, pds_acc) in isolates, compute consensus serovar.

    Returns list of (pathogen, pds_acc, consensus) tuples.
    Does not write to the database — caller decides when/how to persist.
    """
    rows = conn.execute("""
        SELECT pathogen, pds_acc, serovar
        FROM isolates
        WHERE pds_acc IS NOT NULL AND pds_acc != ''
        ORDER BY pathogen, pds_acc
    """).fetchall()

    by_cluster: dict[tuple[str, str], list[str | None]] = {}
    for r in rows:
        key = (r["pathogen"], r["pds_acc"])
        by_cluster.setdefault(key, []).append(r["serovar"])

    out: list[tuple[str, str, SerovarConsensus]] = []
    for (pathogen, pds_acc), serovars in by_cluster.items():
        consensus = compute_consensus_serovar(serovars)
        out.append((pathogen, pds_acc, consensus))

    log.info("Computed consensus serovar for %d clusters", len(out))
    return out
