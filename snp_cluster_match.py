#!/usr/bin/env python3
"""
snp_cluster_match.py

Match an uploaded Salmonella assembly to the nearest precomputed SNP cluster by
cgMLST allele distance, against the 1,503 cluster-representative allelic profiles.

Pipeline this fits into:
  1. chewBBACA AlleleCall on the 1,503 reps (done) -> results_alleles.tsv
  2. chewBBACA ExtractCgMLST                        -> cgMLST95.tsv, cgMLSTschema95.txt
  3. build-index : clean + load the reference matrix, attach SNP-cluster labels from sqlite
  4. calibrate   : within- vs between-cluster distance distributions -> choose a cutoff
  5. query       : call a new assembly (--no-inferred --gl), score it against the index

Allele-distance convention (cgMLST Hamming with missing-data masking):
  reference cell: positive int = allele id ; 0 = missing (special class / LNF)
  query cell:     positive int = allele id (EXC; present in schema -> present in a ref)
                  -1           = INF, novel vs schema. The schema was built from the refs,
                                 so a novel allele cannot match any ref -> counts as a
                                 difference wherever the ref has a call. NOT masked.
                  0            = special class / LNF -> masked
  distance(q, r) = count of core loci where BOTH have a call and the ids differ.
  A locus is scored only if both q and r have a call (zeros mask either side).
"""

import argparse
import json
import os
import re
import sqlite3
from collections import Counter, defaultdict
import subprocess
import sys
import tempfile
from glob import glob

import numpy as np

# ---------------------------------------------------------------------------
# Config. Fill DB_* after running `inspect-db`. Adjust ID handling to match how
# your .fna files are named relative to the sqlite key (accession, PDT, etc.).
# ---------------------------------------------------------------------------
DB_TABLE = "isolates"
DB_ID_COL = "asm_acc"
DB_CLUSTER_COL = "pds_acc"
ID_REGEX = r"(GC[AF]_\d+\.\d+)"
MIN_SHARED_FRAC = 0.80   # below this fraction of core loci shared, a call is flagged low-confidence


# ---------------------------------------------------------------------------
# value cleaning
# ---------------------------------------------------------------------------
_INT_RE = re.compile(r"^\d+$")


def clean_ref_value(s):
    """Reference cell -> int. INF-N -> N (these are real schema alleles). Non-numeric -> 0."""
    s = s.strip()
    if s.startswith("INF-"):
        s = s[4:]
    return int(s) if _INT_RE.match(s) else 0


def clean_query_value(s):
    """Query cell -> int. EXC int -> int. INF* -> -1 (novel, differs from all refs). Else -> 0."""
    s = s.strip()
    if s.startswith("INF"):
        return -1
    return int(s) if _INT_RE.match(s) else 0


def _norm_locus(name):
    """Normalize a locus identifier so schema list and TSV headers align."""
    name = name.strip()
    return name[:-6] if name.endswith(".fasta") else name


# ---------------------------------------------------------------------------
# reference matrix
# ---------------------------------------------------------------------------
def load_reference_matrix(matrix_path, paralogs_path=None):
    """Read a chewBBACA ExtractCgMLST profile (cgMLST95.tsv). Returns (genome_ids, locus_names, R)."""
    drop = set()
    if paralogs_path and os.path.exists(paralogs_path):
        with open(paralogs_path) as fh:
            next(fh, None)  # header
            for line in fh:
                tok = line.split("\t")[0].strip()
                if tok:
                    drop.add(_norm_locus(tok))
        if drop:
            print(f"[build-index] dropping {len(drop)} paralogous loci", file=sys.stderr)

    with open(matrix_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        loci_all = [_norm_locus(h) for h in header[1:]]
        keep_idx = [i for i, loc in enumerate(loci_all) if loc not in drop]
        locus_names = [loci_all[i] for i in keep_idx]

        genome_ids, rows = [], []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genome_ids.append(_norm_locus(parts[0]))  # FILE; strip any .fasta if present
            vals = parts[1:]
            rows.append([clean_ref_value(vals[i]) for i in keep_idx])

    R = np.asarray(rows, dtype=np.int32)
    print(f"[build-index] {R.shape[0]} genomes x {R.shape[1]} core loci", file=sys.stderr)
    return genome_ids, locus_names, R


def load_profiles_aligned(path, locus_names):
    """Read a raw chewBBACA results_alleles.tsv (e.g. a member run) and align its columns to
    locus_names. Loci absent from the file are filled with 0. Cleaned like the reference matrix."""
    with open(path) as fh:
        header = [_norm_locus(h) for h in fh.readline().rstrip("\n").split("\t")[1:]]
        pos = {loc: k for k, loc in enumerate(header)}
        col_for = [pos.get(loc) for loc in locus_names]
        gids, rows = [], []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            gids.append(_norm_locus(parts[0]))
            vals = parts[1:]
            rows.append([clean_ref_value(vals[c]) if (c is not None and c < len(vals)) else 0
                         for c in col_for])
    return gids, np.asarray(rows, dtype=np.int32)


def _db_key(stem):
    if ID_REGEX:
        m = re.search(ID_REGEX, stem)
        return m.group(1) if m else stem
    return stem


def attach_clusters(genome_ids, db_path):
    """Map each genome id -> SNP cluster via sqlite. Unmapped genomes get label ''."""
    if not (DB_TABLE and DB_ID_COL and DB_CLUSTER_COL):
        print("[build-index] DB_TABLE/DB_ID_COL/DB_CLUSTER_COL not set; run inspect-db and edit the "
              "config block. Storing empty cluster labels for now.", file=sys.stderr)
        return np.array(["" for _ in genome_ids], dtype=object)

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    lookup = {}
    for key, clust in cur.execute(f"SELECT {DB_ID_COL}, {DB_CLUSTER_COL} FROM {DB_TABLE}"):
        lookup[str(key)] = "" if clust is None else str(clust)
    con.close()

    clusters, miss = [], 0
    for gid in genome_ids:
        c = lookup.get(_db_key(gid), "")
        if not c:
            miss += 1
        clusters.append(c)
    if miss:
        print(f"[build-index] {miss}/{len(genome_ids)} genomes had no cluster mapping", file=sys.stderr)
    return np.array(clusters, dtype=object)


# ---------------------------------------------------------------------------
# distance
# ---------------------------------------------------------------------------
def masked_hamming(q, R):
    """q:(L,) int, R:(N,L) int. Returns (dist:(N,), shared:(N,))."""
    qp = q != 0                # query has a call (incl. -1 novel)
    rp = R != 0                # ref has a call
    both = rp & qp             # broadcast qp across rows
    diff = both & (R != q)
    return diff.sum(axis=1).astype(np.int32), both.sum(axis=1).astype(np.int32)


# ---------------------------------------------------------------------------
# chewBBACA query call + parse
# ---------------------------------------------------------------------------
def run_allelecall(assembly, schema_dir, core_loci, cpu, workdir):
    indir = os.path.join(workdir, "in")
    outdir = os.path.join(workdir, "out")
    os.makedirs(indir, exist_ok=True)
    dst = os.path.join(indir, os.path.basename(assembly))
    if not os.path.exists(dst):
        os.symlink(os.path.abspath(assembly), dst)
    cmd = ["chewBBACA.py", "AlleleCall", "-i", indir, "-g", schema_dir, "-o", outdir,
           "--gl", core_loci, "--no-inferred", "--cpu", str(cpu)]
    subprocess.run(cmd, check=True)
    hits = glob(os.path.join(outdir, "**", "results_alleles.tsv"), recursive=True)
    if not hits:
        raise FileNotFoundError("results_alleles.tsv not produced by AlleleCall")
    return hits[0]


def parse_query_profile(results_alleles, locus_names):
    """Return a query vector aligned to locus_names (0 where a locus is absent from the call)."""
    with open(results_alleles) as fh:
        header = [_norm_locus(h) for h in fh.readline().rstrip("\n").split("\t")[1:]]
        row = fh.readline().rstrip("\n").split("\t")[1:]
    d = {loc: clean_query_value(v) for loc, v in zip(header, row)}
    return np.asarray([d.get(loc, 0) for loc in locus_names], dtype=np.int32)


# ---------------------------------------------------------------------------
# subcommands
# ---------------------------------------------------------------------------
def cmd_inspect_db(a):
    con = sqlite3.connect(a.db)
    cur = con.cursor()
    tables = [r[0] for r in cur.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
    for t in tables:
        cols = [(r[1], r[2]) for r in cur.execute(f"PRAGMA table_info({t})")]
        n = cur.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        print(f"\n=== {t}  ({n} rows) ===")
        for name, typ in cols:
            print(f"    {name:<32} {typ}")
        for row in cur.execute(f"SELECT * FROM {t} LIMIT 2"):
            print("    sample:", row)
    con.close()


def cmd_build_index(a):
    genome_ids, locus_names, R = load_reference_matrix(a.matrix, a.paralogs)
    if a.members:
        mids, Rm = load_profiles_aligned(a.members, locus_names)
        print(f"[build-index] adding {Rm.shape[0]} member profiles (nearest-member index)",
              file=sys.stderr)
        genome_ids = genome_ids + mids
        R = np.vstack([R, Rm])
    clusters = attach_clusters(genome_ids, a.db) if a.db else np.array(["" for _ in genome_ids], dtype=object)
    np.savez_compressed(a.out, R=R,
                        genome_ids=np.array(genome_ids, dtype=object),
                        locus_names=np.array(locus_names, dtype=object),
                        clusters=clusters)
    print(f"[build-index] wrote {a.out}  ({R.shape[0]} genomes)")


def _load_index(path):
    z = np.load(path, allow_pickle=True)
    return (list(z["genome_ids"]), list(z["locus_names"]),
            z["R"].astype(np.int32), z["clusters"])


def cmd_calibrate(a):
    genome_ids, locus_names, R, clusters = _load_index(a.index)
    N, L = R.shape
    min_shared = int(MIN_SHARED_FRAC * L)
    assigned = np.array([bool(c) for c in clusters])
    clab = np.array([str(c) for c in clusters], dtype=object)
    sizes = Counter(clab[i] for i in range(N) if assigned[i])
    max_size = max(sizes.values()) if sizes else 0

    within, nn_other = [], []          # within-cluster pair dists ; nearest other-cluster dist per genome
    for i in range(N):
        if not assigned[i]:
            continue
        dist, shared = masked_hamming(R[i], R)
        ok = shared >= min_shared
        ok[i] = False
        same = assigned & (clab == clab[i]) & ok
        diff = assigned & (clab != clab[i]) & ok
        for j in np.where(same)[0]:
            if j > i:
                within.append(int(dist[j]))
        if diff.any():
            nn_other.append(int(dist[diff].min()))

    within = np.array(within)
    nn = np.array(nn_other)
    print(f"genomes assigned to a cluster: {int(assigned.sum())}/{N}")
    print(f"distinct clusters: {len(sizes)}   max reps in any one cluster: {max_size}")
    print(f"core loci: {L}   min shared loci for a counted pair: {min_shared}\n")

    if within.size:
        print(f"within-cluster pairs: {within.size}")
        print(f"  dist median/95th/99th/max: "
              f"{np.median(within):.0f} / {np.percentile(within,95):.0f} / "
              f"{np.percentile(within,99):.0f} / {within.max():.0f}")
        print(f"\nsuggested same-cluster cutoff (99th pct within-cluster): {int(np.percentile(within,99))}")
        print("Sanity check against known NCBI cluster definitions; allele distance tracks but is "
              "not identical to SNP-tree distance.")
    else:
        print("within-cluster pairs: 0  (one representative per cluster -> no within-cluster signal)")

    if nn.size:
        print(f"\ninter-cluster separation (nearest other-cluster distance per genome):")
        print(f"  min/1st/5th/median: {nn.min():.0f} / {np.percentile(nn,1):.0f} / "
              f"{np.percentile(nn,5):.0f} / {np.median(nn):.0f}")

    if not within.size:
        print("\nWith one rep per cluster the cutoff must come from a holdout, not this set: pull "
              "non-representative isolates that have an asm_acc and a pds_acc among your 1,503, call "
              "them, and measure distance to their own cluster's rep versus the nearest other rep. "
              "The inter-cluster floor above is the ceiling a real cutoff must sit well below. Treat "
              "'nearest rep' as nearest lineage, not proof of cluster identity: tightest for clonal "
              "clusters, loosest for diverse ones.")


def enrich_cluster(db_path, pds):
    """Pull display metadata for a called cluster from cluster_typing + cluster_summary."""
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    out = {"pds": pds}
    row = cur.execute("SELECT consensus_serovar, mlst_st FROM cluster_typing WHERE pds_acc=?",
                      (pds,)).fetchone()
    if row:
        out["serovar"], out["mlst_st"] = row
    row = cur.execute(
        "SELECT n_total,n_human,n_nonhuman,n_food,n_animal,n_environment,earliest_collection_date,"
        "latest_collection_date,temporal_span_days,countries_json,signals_json "
        "FROM cluster_summary WHERE pds_acc=?", (pds,)).fetchone()
    con.close()
    if row:
        (out["n_total"], out["n_human"], out["n_nonhuman"], out["n_food"], out["n_animal"],
         out["n_env"], out["earliest"], out["latest"], out["span_days"], cj, sj) = row
        try:
            out["top_countries"] = [(d.get("country"), d.get("n")) for d in (json.loads(cj) or [])[:3]]
        except Exception:
            out["top_countries"] = []
        try:
            sig = json.loads(sj) if sj else {}
            out["signals"] = [k for k, v in sig.items() if isinstance(v, dict) and v.get("flag")]
        except Exception:
            out["signals"] = []
    return out


def cmd_query(a):
    genome_ids, locus_names, R, clusters = _load_index(a.index)
    with tempfile.TemporaryDirectory() as wd:
        res = run_allelecall(a.assembly, a.schema, a.core_loci, a.cpu, wd)
        q = parse_query_profile(res, locus_names)
    dist, shared = masked_hamming(q, R)
    L = R.shape[1]
    order = np.lexsort((-shared, dist))  # nearest first, ties broken by more shared loci
    print(f"query: {os.path.basename(a.assembly)}   core loci: {L}\n")
    print(f"{'rank':<5}{'genome':<30}{'cluster':<16}{'dist':>6}{'shared':>8}{'frac':>7}")
    for rank, idx in enumerate(order[:a.topk], 1):
        frac = shared[idx] / L
        flag = "" if frac >= MIN_SHARED_FRAC else "  (low overlap)"
        print(f"{rank:<5}{str(genome_ids[idx]):<30}{str(clusters[idx]) or '-':<16}"
              f"{dist[idx]:>6}{shared[idx]:>8}{frac:>7.2f}{flag}")

    best = order[0]
    best_clu = str(clusters[best])
    # nearest member of a *different* assigned cluster -> the competitor
    second = next((idx for idx in order if str(clusters[idx]) and str(clusters[idx]) != best_clu), None)

    print(f"\nnearest cluster: {best_clu or '(unassigned)'}  dist={dist[best]}  "
          f"(nearest member {genome_ids[best]})")
    if second is not None:
        margin = int(dist[second] - dist[best])
        print(f"next cluster:    {clusters[second]}  dist={dist[second]}  margin={margin}")
        if margin < a.min_margin:
            print(f"[AMBIGUOUS] top two clusters within {margin} alleles (< {a.min_margin}); "
                  f"not a confident single-cluster call")
    if a.threshold is not None and dist[best] > a.threshold:
        print(f"[NO CLOSE MATCH] nearest is {dist[best]} alleles away (> candidate cutoff {a.threshold}); "
              f"may be an unrepresented lineage")
    if shared[best] / L < MIN_SHARED_FRAC:
        print("[LOW OVERLAP] few shared loci; check assembly quality/species before trusting the call")

    if getattr(a, "enrich", False) and a.db and best_clu:
        info = enrich_cluster(a.db, best_clu)
        print(f"\n--- {best_clu} ---")
        print(f"serovar: {info.get('serovar') or 'n/a'}   MLST ST: {info.get('mlst_st') or 'n/a'}")
        if "n_total" in info:
            print(f"isolates: {info['n_total']} (human {info['n_human']}, non-human {info['n_nonhuman']}; "
                  f"food {info['n_food']}, animal {info['n_animal']}, env {info['n_env']})")
            print(f"span: {info['earliest']} to {info['latest']} ({info['span_days']} days)")
            tc = [f"{c} ({n})" for c, n in info.get("top_countries", []) if c]
            if tc:
                print("top countries: " + ", ".join(tc))
            if info.get("signals"):
                print("signals: " + ", ".join(info["signals"]))


def cmd_densify_manifest(a):
    """List additional assemblies (with known cluster) to add to the index, capping per cluster
    and prioritising clusters that are currently thin. Writes a manifest and an accession list."""
    genome_ids, locus_names, R, clusters = _load_index(a.index)
    idx_count = Counter(str(c) for c in clusters if c)
    have_keys = {_db_key(str(g)) for g in genome_ids}

    con = sqlite3.connect(a.db)
    cur = con.cursor()
    pool = defaultdict(list)
    for pds, asm in cur.execute(
            f"SELECT {DB_CLUSTER_COL}, {DB_ID_COL} FROM {DB_TABLE} WHERE {DB_ID_COL} IS NOT NULL"):
        pds = str(pds)
        if pds in idx_count and asm not in have_keys:
            pool[pds].append(asm)
    con.close()

    manifest = []
    # bring each cluster up to --per-cluster total members, thinnest first
    for pds in sorted(idx_count, key=lambda c: idx_count[c]):
        need = a.per_cluster - idx_count[pds]
        if need <= 0:
            continue
        seen = set()
        for asm in pool.get(pds, []):
            if asm in seen:
                continue
            seen.add(asm)
            manifest.append((asm, pds))
            if len(seen) >= need:
                break

    with open(a.out_manifest, "w") as fh:
        fh.write("asm_acc\tpds_acc\n")
        for asm, pds in manifest:
            fh.write(f"{asm}\t{pds}\n")
    with open(a.out_accessions, "w") as fh:
        fh.write("\n".join(asm for asm, _ in manifest) + ("\n" if manifest else ""))

    enriched = len({pds for _, pds in manifest})
    print(f"clusters in index: {len(idx_count)}   target members/cluster: {a.per_cluster}")
    print(f"additional assemblies to fetch: {len(manifest)}  across {enriched} clusters")
    print(f"wrote {a.out_manifest} and {a.out_accessions}")
    print("\nNext: download these (e.g. `datasets download genome accession --inputfile "
          f"{a.out_accessions}`), unzip the .fna into a folder, then run AlleleCall in default mode "
          "(schema-updating) restricted to the core loci:\n"
          "  chewBBACA.py AlleleCall -i members_dir/ -g <working_schema> \\\n"
          "    -o members_out/ --gl cgMLSTschema95.txt --cpu 14\n"
          "Then rebuild: build-index --matrix cgMLST95.tsv --paralogs ... --members "
          "members_out/results_alleles.tsv --db ... --out cluster_index.npz")


def cmd_validate(a):
    """Leave-one-out nearest-neighbour accuracy. Each genome with a same-cluster sibling is scored
    against all others; a call is correct if the nearest (top-1) or any top-k neighbour shares its
    true cluster. This is the empirical accuracy of nearest-member assignment on the index."""
    genome_ids, locus_names, R, clusters = _load_index(a.index)
    N, L = R.shape
    min_shared = int(MIN_SHARED_FRAC * L)
    clab = np.array([str(c) for c in clusters], dtype=object)
    assigned = clab != ""
    sizes = Counter(clab[i] for i in range(N) if assigned[i])
    has_sib = np.array([assigned[i] and sizes[clab[i]] >= 2 for i in range(N)])

    top1 = topk = evaluated = 0
    k = a.topk
    report = getattr(a, "report_misses", False)
    min_margin = getattr(a, "min_margin", 10)
    misses = []
    for i in range(N):
        if not has_sib[i]:
            continue
        dist, shared = masked_hamming(R[i], R)
        ok = (shared >= min_shared) & assigned
        ok[i] = False
        cand = np.where(ok)[0]
        if cand.size == 0:
            continue
        evaluated += 1
        ranked = cand[np.lexsort((-shared[cand], dist[cand]))]
        rclab = clab[ranked]
        true = clab[i]
        if rclab[0] == true:
            top1 += 1
        if true in rclab[:k]:
            topk += 1

        if report and rclab[0] != true:
            pred = rclab[0]
            dist_pred = int(dist[ranked[0]])
            t_hits = np.where(rclab == true)[0]               # where the true cluster first appears
            dist_true = int(dist[ranked[t_hits[0]]]) if t_hits.size else None
            sec = next((j for j in range(ranked.size) if rclab[j] != pred), None)  # 2nd distinct cluster
            second_clu = rclab[sec] if sec is not None else None
            margin_sec = int(dist[ranked[sec]] - dist_pred) if sec is not None else None
            misses.append({
                "genome": str(genome_ids[i]), "true_cluster": true, "pred_cluster": pred,
                "dist_pred": dist_pred, "dist_true": dist_true,
                "gap": (dist_true - dist_pred) if dist_true is not None else None,
                "second_cluster": second_clu, "margin_to_second": margin_sec,
                "true_is_runner_up": (second_clu == true),
                "flagged": (margin_sec is not None and margin_sec < min_margin),
            })

    print(f"validatable genomes (cluster size >= 2): {int(has_sib.sum())}")
    print(f"evaluated (enough shared loci):          {evaluated}")
    if evaluated:
        print(f"top-1 nearest-member accuracy: {top1}/{evaluated} = {top1/evaluated:.1%}")
        print(f"top-{k} recall:                 {topk}/{evaluated} = {topk/evaluated:.1%}")

    if report:
        n = len(misses)
        flagged = sum(m["flagged"] for m in misses)
        runner_up = sum(m["true_is_runner_up"] for m in misses)
        print(f"\n--- miss characterization (--min-margin {min_margin}) ---")
        print(f"misses (top-1 wrong): {n}")
        if n:
            print(f"  caught by ambiguity flag (next cluster within {min_margin}): "
                  f"{flagged}/{n} = {flagged/n:.1%}")
            print(f"  CONFIDENT WRONG (next cluster >= {min_margin} away):        "
                  f"{n-flagged}/{n} = {(n-flagged)/n:.1%}")
            print(f"  true cluster is the immediate runner-up:                    "
                  f"{runner_up}/{n} = {runner_up/n:.1%}")
            gaps = np.array([m["gap"] for m in misses if m["gap"] is not None])
            if gaps.size:
                print(f"  gap, true-cluster dist minus wrong dist, median/95th/max: "
                      f"{np.median(gaps):.0f} / {np.percentile(gaps,95):.0f} / {gaps.max():.0f}")
            cg = np.array([m["gap"] for m in misses if not m["flagged"] and m["gap"] is not None])
            if cg.size:
                print(f"  among confident-wrong, gap median/max: {np.median(cg):.0f} / {cg.max():.0f}")
        out = getattr(a, "miss_out", None)
        if out and misses:
            import csv
            cols = list(misses[0].keys())
            with open(out, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
                w.writeheader()
                w.writerows(misses)
            print(f"  wrote per-miss detail to {out}")
        print("\nFlagged misses are honest uncertainty: the tool would mark them AMBIGUOUS rather than "
              "call them confidently. Confident-wrong is the set that matters; a large gap there means "
              "the genome genuinely sits closer to another cluster on the core than to its own, the "
              "irreducible single-linkage mismatch, not a tunable error.")
        return

    print("\nRun before and after densifying. Top-1 is how often the single nearest member is the "
          "right cluster; top-k recall is how often the right cluster appears at all in the top hits. "
          "Densification should move both up, top-k first.")


def main():
    p = argparse.ArgumentParser(description="Match an assembly to the nearest SNP cluster by cgMLST allele distance")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("inspect-db", help="print sqlite tables/columns to wire the cluster mapping")
    s.add_argument("--db", required=True)
    s.set_defaults(func=cmd_inspect_db)

    s = sub.add_parser("build-index", help="load cgMLST matrix + cluster labels into a .npz index")
    s.add_argument("--matrix", required=True, help="ExtractCgMLST cgMLST95.tsv")
    s.add_argument("--members", help="results_alleles.tsv from a member AlleleCall run (densification)")
    s.add_argument("--db", help="pathogen-watch.sqlite (optional until DB config is filled)")
    s.add_argument("--paralogs", help="paralogous_counts.tsv to exclude its loci")
    s.add_argument("--out", default="cluster_index.npz")
    s.set_defaults(func=cmd_build_index)

    s = sub.add_parser("calibrate", help="within- vs between-cluster distance distributions")
    s.add_argument("--index", required=True)
    s.set_defaults(func=cmd_calibrate)

    s = sub.add_parser("validate", help="leave-one-out nearest-member top-1/top-k accuracy")
    s.add_argument("--index", required=True)
    s.add_argument("--topk", type=int, default=5)
    s.add_argument("--report-misses", dest="report_misses", action="store_true",
                   help="characterize top-1 misses: caught by ambiguity flag vs confident-wrong")
    s.add_argument("--min-margin", dest="min_margin", type=int, default=10,
                   help="margin under which a miss would be flagged AMBIGUOUS (match query default)")
    s.add_argument("--miss-out", dest="miss_out", help="write per-miss detail TSV to this path")
    s.set_defaults(func=cmd_validate)

    s = sub.add_parser("densify-manifest", help="list extra assemblies to fetch per cluster")
    s.add_argument("--index", required=True)
    s.add_argument("--db", required=True)
    s.add_argument("--per-cluster", dest="per_cluster", type=int, default=5,
                   help="target total members per cluster")
    s.add_argument("--out-manifest", dest="out_manifest", default="densify_manifest.tsv")
    s.add_argument("--out-accessions", dest="out_accessions", default="densify_accessions.txt")
    s.set_defaults(func=cmd_densify_manifest)

    s = sub.add_parser("query", help="score one assembly against the index")
    s.add_argument("--index", required=True)
    s.add_argument("--assembly", required=True, help="FASTA (.fna)")
    s.add_argument("--schema", required=True, help="schema dir (working schema, kept current)")
    s.add_argument("--core-loci", dest="core_loci", required=True, help="cgMLSTschema95.txt")
    s.add_argument("--cpu", type=int, default=8)
    s.add_argument("--topk", type=int, default=10)
    s.add_argument("--min-margin", dest="min_margin", type=int, default=10,
                   help="flag AMBIGUOUS if top two clusters are within this many alleles")
    s.add_argument("--threshold", type=int, default=None, help="candidate cutoff for 'no close match'")
    s.add_argument("--enrich", action="store_true", help="print serovar/host/geography/signals for the call")
    s.add_argument("--db", help="pathogen-watch.sqlite, required with --enrich")
    s.set_defaults(func=cmd_query)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
