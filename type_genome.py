#!/usr/bin/env python3
"""
type_genome.py

Run the full SalmonellaWatch-Ag typing panel on one assembly:
  - SNP cluster call (snp_cluster_match.py query --enrich)
  - Serovar (SISTR)
  - 7-gene Achtman MLST (mlst, senterica scheme)
  - Plasmid replicons (abricate, PlasmidFinder database)
  - pESI backbone markers (pesi_markers.py)
  - Virulence & AMR (AMRFinderPlus, --plus --organism Salmonella)

Emits a single JSON record per genome. Each tool's contribution is independent: if
a tool isn't installed or fails, that section reports the error and the rest still
runs. Designed to be called from a shell wrapper that resolves paths.

Honest scope: this is the in silico surveillance profile of a draft assembly. Every
section depends on assembly quality; the pESI section in particular is a marker-
supported inference, not a confirmed plasmid, and any positive call there requires
manual review.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from glob import glob


def _which(name):
    return shutil.which(name) is not None


def _run(cmd, cwd=None, check=True):
    """Run a command, return (returncode, stdout, stderr)."""
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if check and p.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {p.stderr.strip()[:500]}")
    return p.returncode, p.stdout, p.stderr


def section_cluster(assembly, args):
    if not args.index or not args.schema or not args.core_loci:
        return {"status": "skipped", "reason": "--index/--schema/--core-loci not provided"}
    cmd = [sys.executable, os.path.join(args.script_dir, "snp_cluster_match.py"), "query",
           "--index", args.index, "--assembly", assembly, "--schema", args.schema,
           "--core-loci", args.core_loci, "--cpu", str(args.cpu)]
    if args.db:
        cmd += ["--enrich", "--db", args.db]
    try:
        _, out, _ = _run(cmd)
    except RuntimeError as e:
        return {"status": "error", "error": str(e)}
    # parse the verdict and enrich block from stdout
    info = {"status": "ok", "raw": out.strip()}
    m = re.search(r"nearest cluster: (\S+)\s+dist=(\d+)", out)
    if m:
        info["pds"] = m.group(1); info["distance"] = int(m.group(2))
    m = re.search(r"next cluster:\s+(\S+)\s+dist=(\d+)\s+margin=(\d+)", out)
    if m:
        info["next_pds"] = m.group(1); info["next_distance"] = int(m.group(2))
        info["margin"] = int(m.group(3))
    info["ambiguous"] = "[AMBIGUOUS]" in out
    info["no_close_match"] = "[NO CLOSE MATCH]" in out
    info["low_overlap"] = "[LOW OVERLAP]" in out
    m = re.search(r"serovar:\s+(.+?)\s{3}MLST ST: (.+?)$", out, re.M)
    if m:
        info["enrich"] = {"serovar": m.group(1).strip(), "mlst_st_cluster": m.group(2).strip()}
    return info


def section_sistr(assembly):
    if not _which("sistr"):
        return {"status": "skipped", "reason": "sistr not on PATH"}
    with tempfile.TemporaryDirectory() as wd:
        out_json = os.path.join(wd, "out.json")
        try:
            _run(["sistr", "-f", "json", "-o", out_json, "--no-cgmlst", assembly])
            with open(out_json) as fh:
                data = json.load(fh)
            d = data[0] if isinstance(data, list) and data else data
            return {"status": "ok", "serovar": d.get("serovar"),
                    "serogroup": d.get("serogroup"),
                    "h1": d.get("h1"), "h2": d.get("h2"),
                    "antigenic_formula": d.get("serovar_antigen"),
                    "qc_status": d.get("qc_status"), "qc_messages": d.get("qc_messages")}
        except Exception as e:
            return {"status": "error", "error": str(e)}


def section_mlst(assembly):
    if not _which("mlst"):
        return {"status": "skipped", "reason": "mlst not on PATH"}
    try:
        _, out, _ = _run(["mlst", "--scheme", "senterica_achtman_2", "--quiet", assembly])
    except RuntimeError:
        # fall back to auto-detect if Achtman scheme name differs in user's mlst install
        try:
            _, out, _ = _run(["mlst", "--quiet", assembly])
        except RuntimeError as e:
            return {"status": "error", "error": str(e)}
    parts = out.strip().split("\t")
    if len(parts) < 3:
        return {"status": "error", "error": f"unexpected mlst output: {out[:200]}"}
    return {"status": "ok", "scheme": parts[1], "st": parts[2], "alleles": parts[3:]}


def section_plasmidfinder(assembly, cpu):
    if not _which("abricate"):
        return {"status": "skipped", "reason": "abricate not on PATH"}
    try:
        _, out, _ = _run(["abricate", "--db", "plasmidfinder", "--quiet",
                          "--threads", str(cpu), assembly])
    except RuntimeError as e:
        return {"status": "error", "error": str(e)}
    lines = [l for l in out.splitlines() if l and not l.startswith("#")]
    hits = []
    for line in lines:
        v = line.split("\t")
        if len(v) >= 11:
            hits.append({"replicon": v[5], "pct_cov": float(v[9]), "pct_id": float(v[10])})
    return {"status": "ok", "n_replicons": len(hits), "replicons": hits}


def section_pesi(assembly, script_dir, cpu):
    markers_fa = os.path.join(script_dir, "pesi_markers.fasta")
    if not os.path.exists(markers_fa):
        return {"status": "skipped", "reason": f"pesi_markers.fasta missing; run pesi_markers.py fetch-refs"}
    try:
        _, out, _ = _run([sys.executable, os.path.join(script_dir, "pesi_markers.py"),
                          "run", "--assembly", assembly, "--markers", markers_fa, "--cpu", str(cpu)])
    except RuntimeError as e:
        return {"status": "error", "error": str(e)}
    # Match only known marker names to avoid catching the table header ("marker") or footer text
    known = {"rdA", "pilL", "sogS", "trbA", "ipf", "ipr2", "IncFIB_pN55391"}
    present, all_seen = [], []
    for line in out.splitlines():
        m = re.match(r"^(\S+)\s+(present|below|absent)\b", line)
        if not m or m.group(1) not in known:
            continue
        all_seen.append(m.group(1))
        if m.group(2) == "present":
            present.append(m.group(1))
    vmatch = re.search(r"verdict: (.+)$", out, re.M)
    n_in_fasta = 0
    try:
        with open(markers_fa) as fh:
            n_in_fasta = sum(1 for l in fh if l.startswith(">"))
    except OSError:
        pass
    return {"status": "ok", "markers_present": present, "n_present": len(present),
            "n_markers_in_panel": n_in_fasta, "n_total_possible": len(known),
            "verdict": vmatch.group(1) if vmatch else None,
            "caveat": "marker-supported inference on draft assembly; manual review required"}


def section_amrfinderplus(assembly, cpu):
    if not _which("amrfinder"):
        return {"status": "skipped", "reason": "amrfinder not on PATH"}
    try:
        _, out, _ = _run(["amrfinder", "--nucleotide", assembly, "--organism", "Salmonella",
                          "--plus", "--threads", str(cpu)])
    except RuntimeError as e:
        return {"status": "error", "error": str(e)}
    lines = out.splitlines()
    if not lines:
        return {"status": "ok", "amr": [], "virulence": [], "stress": [], "point_mutations": []}
    header = lines[0].split("\t")
    cols = {h: i for i, h in enumerate(header)}
    def col(row, name):
        i = cols.get(name)
        return row[i] if i is not None and i < len(row) else None
    buckets = {"AMR": [], "VIRULENCE": [], "STRESS": [], "POINT": []}
    for line in lines[1:]:
        if not line.strip():
            continue
        r = line.split("\t")
        typ = col(r, "Type") or ""
        rec = {"gene": col(r, "Element symbol"), "name": col(r, "Element name"),
               "subtype": col(r, "Subtype"), "method": col(r, "Method"),
               "pct_cov": _safe_float(col(r, "% Coverage of reference")),
               "pct_id": _safe_float(col(r, "% Identity to reference")),
               "contig": col(r, "Contig id")}
        sub = (col(r, "Subtype") or "").upper()
        if "POINT" in sub:
            buckets["POINT"].append(rec)
        elif "VIRULENCE" in typ.upper() or "VIRULENCE" in sub:
            buckets["VIRULENCE"].append(rec)
        elif "STRESS" in typ.upper() or "STRESS" in sub:
            buckets["STRESS"].append(rec)
        else:
            buckets["AMR"].append(rec)
    return {"status": "ok", "amr": buckets["AMR"], "virulence": buckets["VIRULENCE"],
            "stress": buckets["STRESS"], "point_mutations": buckets["POINT"]}


def _safe_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description="Run full typing panel on one assembly")
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--out", help="write JSON to this path (default: stdout)")
    ap.add_argument("--script-dir", dest="script_dir", default=os.path.dirname(os.path.abspath(__file__)),
                    help="where snp_cluster_match.py and pesi_markers.py live")
    # cluster-call inputs (optional; section is skipped if any missing)
    ap.add_argument("--index"); ap.add_argument("--schema"); ap.add_argument("--core-loci", dest="core_loci")
    ap.add_argument("--db", help="pathogen-watch.sqlite for cluster --enrich")
    # which sections to skip
    ap.add_argument("--skip", nargs="*", default=[],
                    choices=["cluster", "sistr", "mlst", "plasmid", "pesi", "amrfinder"])
    args = ap.parse_args()

    if not os.path.exists(args.assembly):
        sys.exit(f"assembly not found: {args.assembly}")

    sections = {}
    if "cluster" not in args.skip:    sections["cluster"]    = section_cluster(args.assembly, args)
    if "sistr" not in args.skip:      sections["sistr"]      = section_sistr(args.assembly)
    if "mlst" not in args.skip:       sections["mlst"]       = section_mlst(args.assembly)
    if "plasmid" not in args.skip:    sections["plasmid"]    = section_plasmidfinder(args.assembly, args.cpu)
    if "pesi" not in args.skip:       sections["pesi"]       = section_pesi(args.assembly, args.script_dir, args.cpu)
    if "amrfinder" not in args.skip:  sections["amrfinder"]  = section_amrfinderplus(args.assembly, args.cpu)

    record = {"assembly": os.path.basename(args.assembly), "sections": sections}
    payload = json.dumps(record, indent=2)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(payload + "\n")
    else:
        print(payload)


if __name__ == "__main__":
    main()
