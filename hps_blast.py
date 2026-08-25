#!/usr/bin/env python3
"""
hps_blast.py

Detect the 8 Harhay HPS virulence genes by BLASTing their SL1344 reference CDS
against an assembly. More sensitive than in silico PCR for gene carriage: tolerant
of point mutations, and union HSP coverage means a gene split across contigs still
scores as present.

Two modes:
  fetch-refs : download SL1344 (GCF_000210855.2) and extract the 8 CDS by locus tag
  run        : blastn the reference CDS against an assembly and report presence

References are the SL1344_RS locus tags from Harhay et al. 2025 Table 2. SL1344 is a
Typhimurium strain, so for divergent serovars set thresholds with that in mind: a
strict identity cutoff can miss true homologs that differ from the Typhimurium allele.
invA is highly conserved (universal control); the plasmid and effector genes vary more.
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict
from glob import glob

# gene -> SL1344 RefSeq locus tag suffix (Harhay 2025 Table 2)
LOCUS_MAP = {
    "sseK2": "RS10980", "sseK3": "RS10000", "avrA": "RS14835", "lpfB": "RS18760",
    "spvD": "RS24010", "sspH2": "RS11515", "gtgA": "RS05020", "invA": "RS14990",
}
SL1344_ACC = "GCF_000210855.2"
_LT = re.compile(r"\[locus_tag=([^\]]+)\]")


def extract_cds(cds_path, locus_map):
    """Pull CDS whose locus tag suffix matches the panel. Returns {gene: sequence}."""
    want = {suf: gene for gene, suf in locus_map.items()}
    out, cur_gene, cur = {}, None, []
    with open(cds_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_gene:
                    out[cur_gene] = "".join(cur)
                cur_gene, cur = None, []
                m = _LT.search(line)
                if m:
                    suffix = m.group(1).split("_")[-1]
                    if suffix in want:
                        cur_gene = want[suffix]
            elif cur_gene:
                cur.append(line.strip())
    if cur_gene:
        out[cur_gene] = "".join(cur)
    return out


def summarize_blast(rows, min_id, min_cov, noise_id=70.0):
    """Best single locus per gene. HSPs below noise_id are dropped, the rest are grouped
    by contig and unioned within each contig, and the contig with the highest coverage x
    identity score is reported. This finds the true ortholog locus instead of summing
    scattered partial/paralog hits from across the genome."""
    by_gene = defaultdict(list)
    for r in rows:
        if r["pident"] >= noise_id:
            by_gene[r["qseqid"]].append(r)
    out = []
    for gene, hsps in by_gene.items():
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
            covered = sum(e - s + 1 for s, e in merged)
            cov = 100.0 * covered / qlen if qlen else 0.0
            tot = sum(h["length"] for h in chs)
            wid = sum(h["pident"] * h["length"] for h in chs) / tot if tot else 0.0
            score = cov * wid                       # favour high coverage AND high identity
            if best is None or score > best["score"]:
                best = {"contig": contig, "cov": cov, "id": wid, "nh": len(chs), "score": score}
        out.append({"gene": gene, "pct_id": round(best["id"], 1), "pct_cov": round(best["cov"], 1),
                    "n_hsps": best["nh"], "contig": best["contig"],
                    "present": best["cov"] >= min_cov and best["id"] >= min_id})
    return out


def run_blastn(genes_fasta, assembly, cpu):
    cols = ["qseqid", "sseqid", "pident", "length", "qlen", "qstart", "qend",
            "sstart", "send", "evalue", "bitscore"]
    cmd = ["blastn", "-task", "dc-megablast", "-query", genes_fasta, "-subject", assembly,
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


def cmd_fetch_refs(a):
    workdir = "sl1344_ref"
    os.makedirs(workdir, exist_ok=True)
    zf = os.path.join(workdir, "sl1344.zip")
    subprocess.run(["datasets", "download", "genome", "accession", SL1344_ACC,
                    "--include", "cds", "--filename", zf], check=True)
    subprocess.run(["unzip", "-o", "-q", zf, "-d", workdir], check=True)
    hits = glob(os.path.join(workdir, "**", "cds_from_genomic.fna"), recursive=True)
    if not hits:
        sys.exit("could not find cds_from_genomic.fna in the download")
    genes = extract_cds(hits[0], LOCUS_MAP)
    with open(a.out, "w") as fh:
        for gene, seq in genes.items():
            fh.write(f">{gene}\n{seq}\n")
    found = set(genes)
    missing = set(LOCUS_MAP) - found
    print(f"wrote {len(found)} reference CDS to {a.out}: {', '.join(sorted(found))}")
    if missing:
        print(f"WARNING: not found by locus tag (check annotation version): {', '.join(sorted(missing))}")


def cmd_run(a):
    rows = run_blastn(a.genes, a.assembly, a.cpu)
    res = {r["gene"]: r for r in summarize_blast(rows, a.min_id, a.min_cov)}
    print(f"thresholds: identity >= {a.min_id}%, coverage >= {a.min_cov}%\n")
    print(f"{'gene':<10}{'present':<9}{'%id':>6}{'%cov':>7}{'hsps':>6}  best contig")
    present = 0
    for gene in LOCUS_MAP:                      # stable panel order, including absent genes
        r = res.get(gene)
        if r is None:
            print(f"{gene:<10}{'absent':<9}{'-':>6}{'-':>7}{'-':>6}  (no hit)")
            continue
        present += int(r["present"])
        tag = "present" if r["present"] else "below"
        print(f"{gene:<10}{tag:<9}{r['pct_id']:>6}{r['pct_cov']:>7}{r['n_hsps']:>6}  {r['contig']}")
    hits = [g for g in LOCUS_MAP if res.get(g, {}).get("present")]
    print(f"\ndetected {present}/{len(LOCUS_MAP)}: {', '.join(hits) if hits else 'none'}")
    print("invA is the universal control; absence suggests a poor or non-Salmonella assembly.")


def main():
    ap = argparse.ArgumentParser(description="BLAST-based HPS virulence gene detection")
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("fetch-refs", help="download SL1344 and extract the 8 CDS by locus tag")
    s.add_argument("--out", default="hps_genes.fasta")
    s.set_defaults(func=cmd_fetch_refs)

    s = sub.add_parser("run", help="blastn reference CDS against an assembly")
    s.add_argument("--assembly", required=True)
    s.add_argument("--genes", default="hps_genes.fasta")
    s.add_argument("--min-id", dest="min_id", type=float, default=90.0)
    s.add_argument("--min-cov", dest="min_cov", type=float, default=90.0)
    s.add_argument("--cpu", type=int, default=4)
    s.set_defaults(func=cmd_run)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
