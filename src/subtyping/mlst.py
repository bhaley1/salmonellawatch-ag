"""Run MLST against an assembly FASTA.

Wraps Torsten Seemann's `mlst` CLI (https://github.com/tseemann/mlst):
a well-established tool that BLASTs an assembly against PubMLST scheme
databases and returns the sequence type.

Usage from the shell:
    $ mlst --scheme lmonocytogenes_2 /path/to/assembly.fna
    /path/to/assembly.fna  lmonocytogenes_2  ST6  abcZ(3)  bglA(1)  ...

Output is tab-delimited: file, scheme, ST, then allele(allele_number) pairs.

We parse that into MLSTResult. If `mlst` is not installed, this module
returns a clear "tool_not_installed" error and the caller decides what
to do (typically: log it, record the error in cluster_typing, move on).
"""

from __future__ import annotations

import json
import logging
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)


# Map our display pathogen names to mlst scheme names.
# The PubMLST scheme catalog uses these identifiers.
PATHOGEN_SCHEME = {
    "Listeria": "lmonocytogenes",
    "Salmonella": "senterica",
    "STEC":     "ecoli_achtman_4",  # Phase 4 will need this
}


@dataclass
class MLSTResult:
    """Outcome of an MLST run against one assembly."""
    scheme: str | None      # e.g. 'lmonocytogenes_2'
    st: str | None          # e.g. 'ST6', 'novel', 'untypeable'
    alleles: dict[str, str] # {gene: allele_number or '~allele' for partial}
    error: str | None       # 'tool_not_installed' | 'no_match' | 'subprocess_error' | None


def is_mlst_installed() -> bool:
    """Whether the `mlst` CLI is on PATH."""
    return shutil.which("mlst") is not None


def run_mlst(fasta_path: Path, scheme: str, timeout: int = 300) -> MLSTResult:
    """Run `mlst` against a FASTA and parse the result.

    Returns an MLSTResult. Never raises — errors are encoded in the result.
    """
    if not is_mlst_installed():
        log.warning("mlst CLI not installed; skipping typing")
        return MLSTResult(
            scheme=scheme, st=None, alleles={}, error="tool_not_installed",
        )

    if not fasta_path.exists():
        return MLSTResult(
            scheme=scheme, st=None, alleles={},
            error=f"fasta_missing:{fasta_path}",
        )

    try:
        proc = subprocess.run(
            ["mlst", "--scheme", scheme, "--quiet", str(fasta_path)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        log.warning("mlst timed out on %s", fasta_path.name)
        return MLSTResult(scheme=scheme, st=None, alleles={}, error="timeout")
    except Exception as e:
        log.warning("mlst subprocess error on %s: %s", fasta_path.name, e)
        return MLSTResult(scheme=scheme, st=None, alleles={}, error=f"subprocess_error:{e}")

    if proc.returncode != 0:
        log.warning("mlst returned %d: stderr=%s", proc.returncode, proc.stderr[:200])
        return MLSTResult(
            scheme=scheme, st=None, alleles={},
            error=f"nonzero_exit:{proc.returncode}",
        )

    return parse_mlst_output(proc.stdout, scheme)


def parse_mlst_output(stdout: str, requested_scheme: str) -> MLSTResult:
    """Parse mlst's tab-separated output.

    Format:  filename<TAB>scheme<TAB>ST<TAB>gene1(allele)<TAB>gene2(allele)...

    ST values:
      - integer N  → ST{N}
      - '-'        → no scheme matched (untypeable)
      - integer with one or more '~' alleles → novel
    """
    lines = [ln.strip() for ln in stdout.splitlines() if ln.strip()]
    if not lines:
        return MLSTResult(
            scheme=requested_scheme, st=None, alleles={}, error="empty_output",
        )

    parts = lines[0].split("\t")
    if len(parts) < 3:
        return MLSTResult(
            scheme=requested_scheme, st=None, alleles={},
            error=f"malformed_output:{lines[0][:120]}",
        )

    _filename, scheme, st_raw = parts[0], parts[1], parts[2]
    allele_calls = parts[3:]

    # Parse allele calls: e.g. "abcZ(3)" or "abcZ(~3)" (partial)
    alleles: dict[str, str] = {}
    for call in allele_calls:
        m = re.match(r"([A-Za-z][A-Za-z0-9_]*)\(([^)]+)\)", call)
        if m:
            gene, allele = m.groups()
            alleles[gene] = allele

    # Normalize the ST string
    if st_raw == "-":
        st = "untypeable"
    elif "~" in st_raw or any("~" in v for v in alleles.values()):
        # Novel allele(s) → strain is novel
        st = f"novel (closest {st_raw.rstrip('~')})"
    else:
        # Either integer or string; prefix with 'ST' for display
        if st_raw.isdigit():
            st = f"ST{st_raw}"
        else:
            st = st_raw

    return MLSTResult(
        scheme=scheme if scheme != "-" else requested_scheme,
        st=st,
        alleles=alleles,
        error=None if st != "untypeable" else "no_match",
    )


def alleles_to_json(alleles: dict[str, str]) -> str:
    return json.dumps(alleles, sort_keys=True)
