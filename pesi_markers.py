#!/usr/bin/env python3
"""
pesi_markers.py

Detect the seven pESI/pESI-like megaplasmid backbone markers from Lee et al. 2024
(PLOS ONE; pN55391 reference): rdA, pilL, SogS, TrbA, ipf, ipr2, IncFIB(pN55391).

Scoring follows Lee et al.'s operational definition: >= 3 markers = pESI-positive;
6-7 markers = high confidence; 4-5 = moderate; <= 2 = unlikely or fragmentary.

CAVEAT: pESI is a ~280-300 kb plasmid. On a short-read draft assembly it
fragments across many contigs, so every in silico call here is a marker-supported
inference, not a confirmed plasmid. ALL positive calls require manual inspection
of the plasmid contigs and, ideally, a long-read or hybrid assembly to confirm.

Modes:
  fetch-refs : pull pN55391 CDS and extract the marker set by gene/product
  run        : blastn the marker fasta against an assembly, score, classify
"""

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from glob import glob

# pN55391 (Salmonella Infantis pESI prototype). NCBI nucleotide accession.
PN55391_NUC = "CP016411.1"
# Marker gene symbols and the synonyms/product hints used to extract them from a
# GenBank/RefSeq CDS dump. The replicon marker is treated specially via repA + IncFIB.
MARKERS = [
    {"name": "rdA",   "aliases": ["ardA"], "product_re": r"antirestriction|ArdA|rdA"},
    {"name": "pilL",  "aliases": ["pilL"], "product_re": r"pilL|conjugal transfer pilus"},
    {"name": "sogS",  "aliases": ["sogS", "SogS"], "product_re": r"SogS|primase"},
    {"name": "trbA",  "aliases": ["trbA", "TrbA"], "product_re": r"TrbA|conjugal transfer"},
    {"name": "ipf",   "aliases": ["ipfA", "ipfB", "ipfC", "ipfD"], "product_re": r"Ipf fimbria"},
    {"name": "ipr2",  "aliases": ["ipr2"], "product_re": r"ipr2"},
    {"name": "IncFIB_pN55391", "aliases": ["repA"], "product_re": r"replication initiat|RepA"},
]

_LT = re.compile(r"\[locus_tag=([^\]]+)\]")
_GENE = re.compile(r"\[gene=([^\]]+)\]")
_PROT = re.compile(r"\[protein=([^\]]+)\]")


def extract_pesi_refs(cds_path):
    """Pull one representative CDS per marker, using gene name first then product text."""
    out = {m["name"]: None for m in MARKERS}
    cur_header, cur_seq = None, []
    records = []
    with open(cds_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_header is not None:
                    records.append((cur_header, "".join(cur_seq)))
                cur_header, cur_seq = line.strip(), []
            else:
                cur_seq.append(line.strip())
    if cur_header is not None:
        records.append((cur_header, "".join(cur_seq)))

    for marker in MARKERS:
        name = marker["name"]
        if out[name]:
            continue
        for header, seq in records:
            g = _GENE.search(header)
            p = _PROT.search(header)
            gene = g.group(1).lower() if g else ""
            prod = p.group(1) if p else ""
            if any(a.lower() == gene for a in marker["aliases"]) or \
               re.search(marker["product_re"], prod, re.IGNORECASE):
                out[name] = seq
                break
    return out


def cmd_fetch_refs(a):
    wd = "pesi_ref"
    os.makedirs(wd, exist_ok=True)
    # The pESI replicon (pN55391) is on a Salmonella Infantis assembly; fetching the
    # nucleotide directly with `datasets` requires the genome accession. Use efetch
    # via entrez-direct (pulled in by the AMRFinderPlus env) for the nucleotide+CDS.
    cds_path = os.path.join(wd, "pn55391_cds.fna")
    if not os.path.exists(cds_path):
        cmd = ["efetch", "-db", "nuccore", "-id", PN55391_NUC, "-format", "fasta_cds_na"]
        with open(cds_path, "w") as out:
            subprocess.run(cmd, check=True, stdout=out)
    refs = extract_pesi_refs(cds_path)
    found = {k: v for k, v in refs.items() if v}
    with open(a.out, "w") as fh:
        for name, seq in found.items():
            fh.write(f">{name}\n{seq}\n")
    print(f"wrote {len(found)}/{len(MARKERS)} pESI marker CDS to {a.out}: {', '.join(found)}")
    missing = [m["name"] for m in MARKERS if m["name"] not in found]
    if missing:
        print(f"NOT FOUND in pN55391 annotation: {', '.join(missing)}", file=sys.stderr)
        print("These need manual sourcing (e.g. from Aviv 2014 pESI annotation).", file=sys.stderr)


def run_blastn(query, assembly, cpu):
    cols = ["qseqid", "sseqid", "pident", "length", "qlen", "qstart", "qend"]
    cmd = ["blastn", "-task", "dc-megablast", "-query", query, "-subject", assembly,
           "-evalue", "1e-10", "-outfmt", "6 " + " ".join(cols), "-num_threads", str(cpu)]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    rows = []
    for line in res.stdout.splitlines():
        if not line.strip():
            continue
        v = line.split("\t")
        rows.append({"qseqid": v[0], "sseqid": v[1], "pident": float(v[2]),
                     "length": int(v[3]), "qlen": int(v[4]),
                     "qstart": int(v[5]), "qend": int(v[6])})
    return rows


def summarize(rows, min_id, min_cov, noise_id=70.0):
    """Best single locus per marker, union coverage within a contig (same logic as the
    HPS detector). Returns one row per marker; markers with no hit are absent."""
    by_marker = defaultdict(list)
    for r in rows:
        if r["pident"] >= noise_id:
            by_marker[r["qseqid"]].append(r)
    out = {}
    for marker, hsps in by_marker.items():
        qlen = hsps[0]["qlen"]
        by_contig = defaultdict(list)
        for h in hsps:
            by_contig[h["sseqid"]].append(h)
        best = None
        for contig, chs in by_contig.items():
            ivs = sorted((min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])) for h in chs)
            merged = []
            for s, e in ivs:
                if merged and s <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            cov = 100.0 * sum(e - s + 1 for s, e in merged) / qlen if qlen else 0.0
            tot = sum(h["length"] for h in chs)
            wid = sum(h["pident"] * h["length"] for h in chs) / tot if tot else 0.0
            score = cov * wid
            if best is None or score > best["score"]:
                best = {"contig": contig, "cov": cov, "id": wid, "score": score}
        out[marker] = {"pct_id": round(best["id"], 1), "pct_cov": round(best["cov"], 1),
                       "contig": best["contig"],
                       "present": best["cov"] >= min_cov and best["id"] >= min_id}
    return out


def classify(n_present):
    """Lee et al. 2024 operational thresholds. >= 3 markers = pESI-positive (99.2%
    PPV in their 8,290-genome cluster). 6-7 = high confidence; 4-5 = moderate."""
    if n_present >= 6:
        return "LIKELY pESI (high confidence)"
    if n_present >= 3:
        return "POSSIBLE pESI (moderate confidence; manual inspection required)"
    if n_present >= 1:
        return "UNLIKELY pESI (partial markers; may reflect fragmentary plasmid or unrelated homologs)"
    return "no pESI markers detected"


def cmd_run(a):
    rows = run_blastn(a.markers, a.assembly, a.cpu)
    res = summarize(rows, a.min_id, a.min_cov)
    print(f"thresholds: identity >= {a.min_id}%, coverage >= {a.min_cov}%\n")
    print(f"{'marker':<20}{'present':<9}{'%id':>6}{'%cov':>7}  contig")
    n_present = 0
    for m in MARKERS:
        name = m["name"]
        r = res.get(name)
        if r is None:
            print(f"{name:<20}{'absent':<9}{'-':>6}{'-':>7}  (no hit)")
            continue
        n_present += int(r["present"])
        tag = "present" if r["present"] else "below"
        print(f"{name:<20}{tag:<9}{r['pct_id']:>6}{r['pct_cov']:>7}  {r['contig']}")
    verdict = classify(n_present)
    print(f"\nmarkers present: {n_present}/{len(MARKERS)}")
    print(f"verdict: {verdict}")
    print("\n* Draft-assembly inference only. Any positive call must be confirmed by")
    print("  manual inspection of the plasmid contigs, ideally with a long-read or")
    print("  hybrid assembly. pESI is ~280-300 kb and routinely fragments across many")
    print("  short-read contigs. Reference plasmid: pN55391 (CP016411.1).")


def main():
    ap = argparse.ArgumentParser(description="pESI backbone marker detection")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("fetch-refs", help="extract the 7 marker CDS from pN55391")
    s.add_argument("--out", default="pesi_markers.fasta")
    s.set_defaults(func=cmd_fetch_refs)
    s = sub.add_parser("run", help="blastn markers vs assembly and score")
    s.add_argument("--assembly", required=True)
    s.add_argument("--markers", default="pesi_markers.fasta")
    s.add_argument("--min-id", dest="min_id", type=float, default=90.0)
    s.add_argument("--min-cov", dest="min_cov", type=float, default=60.0,
                   help="coverage threshold; pESI markers often partially recovered (default 60)")
    s.add_argument("--cpu", type=int, default=4)
    s.set_defaults(func=cmd_run)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
