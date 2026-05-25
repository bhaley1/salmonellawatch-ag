"""SISTR-based serotype prediction for Salmonella.

Runs sistr --no-cgmlst on a representative assembly and returns a
standardized serotype string. The --no-cgmlst flag skips the slow
cgMLST step and uses antigen gene detection only, which is sufficient
for serotype assignment and runs in ~5-10 seconds per assembly.

Within a SNP cluster, all isolates share the same serotype, so we only
need to type the representative assembly once per cluster.
"""

from __future__ import annotations

import logging
import subprocess
import csv
import io
from dataclasses import dataclass
from pathlib import Path

log = logging.getLogger(__name__)


@dataclass
class SISTRResult:
    serovar: str | None
    antigenic_formula: str | None
    serogroup: str | None
    h1: str | None
    h2: str | None
    o_antigen: str | None
    error: str | None = None


# Serotypes of concern from Wheeler et al. 2024 (Front. Microbiol.)
# doi: 10.3389/fmicb.2024.1307563
# Exact SOC list from Wheeler et al. 2024 (Front. Microbiol. doi:10.3389/fmicb.2024.1307563)
# 21 serotypes — do not add or remove without updating the reference
SEROTYPES_OF_CONCERN: set[str] = {
    "Enteritidis",
    "Typhimurium",
    "Heidelberg",
    "Infantis",
    "Newport",
    "Uganda",
    "Braenderup",
    "Muenchen",
    "Montevideo",
    "Javiana",
    "Reading",
    "Dublin",
    "Oranienburg",
    "Potsdam",
    "Thompson",
    "Saintpaul",
    "Hadar",
    "Schwarzengrund",
    "Anatum",
    "Berta",
}

# Antigenic formulas that are SOC even if serovar name differs
SOC_ANTIGENIC_FORMULAS: set[str] = {
    "4:i:-",       # monophasic Typhimurium
    "4,[5]:i:-",
    "4,5:i:-",
    "4,5,12:i:-",
    "1,4,[5],12:i:-",
}


def is_serotype_of_concern(result: SISTRResult) -> bool:
    """Return True if this SISTR result matches a SOC."""
    if result.serovar:
        # Handle pipe-separated predictions like "Typhimurium|Lagos"
        serovars = [s.strip() for s in result.serovar.split("|")]
        for sv in serovars:
            if sv in SEROTYPES_OF_CONCERN:
                return True
    if result.antigenic_formula:
        formula = result.antigenic_formula.strip()
        for soc_formula in SOC_ANTIGENIC_FORMULAS:
            if soc_formula in formula:
                return True
    return False


def clean_serovar(raw: str | None) -> str | None:
    """Return the primary serovar from a pipe-separated SISTR prediction."""
    if not raw:
        return None
    # Take first prediction before pipe
    primary = raw.split("|")[0].strip()
    if not primary or primary in ("-", ""):
        return None
    return primary


def run_sistr(fasta_path: Path, timeout: int = 120) -> SISTRResult:
    """Run sistr --no-cgmlst on a FASTA and return parsed result."""
    import tempfile, os
    tmp_out = tempfile.mktemp(suffix="_sistr")
    try:
        result = subprocess.run(
            ["sistr", "-i", str(fasta_path), "genome",
             "-f", "tab", "-o", tmp_out, "--no-cgmlst"],
            capture_output=True, text=True, timeout=timeout
        )
        tab_file = tmp_out + ".tab"
        stdout = open(tab_file).read() if os.path.exists(tab_file) else ""
        try: os.unlink(tab_file)
        except: pass
    except FileNotFoundError:
        return SISTRResult(
            serovar=None, antigenic_formula=None, serogroup=None,
            h1=None, h2=None, o_antigen=None, error="tool_not_installed"
        )
    except subprocess.TimeoutExpired:
        return SISTRResult(
            serovar=None, antigenic_formula=None, serogroup=None,
            h1=None, h2=None, o_antigen=None, error="timeout"
        )
    except Exception as e:
        return SISTRResult(
            serovar=None, antigenic_formula=None, serogroup=None,
            h1=None, h2=None, o_antigen=None, error=f"error:{e}"
        )

    if result.returncode != 0:
        return SISTRResult(
            serovar=None, antigenic_formula=None, serogroup=None,
            h1=None, h2=None, o_antigen=None,
            error=f"nonzero_exit:{result.returncode}"
        )

    # Parse tab output — header row then data row
    # stdout was read from temp file above
    # Filter out warning lines
    lines = [l for l in stdout.splitlines() if not l.startswith("/") and l.strip()]
    if len(lines) < 2:
        return SISTRResult(
            serovar=None, antigenic_formula=None, serogroup=None,
            h1=None, h2=None, o_antigen=None, error="no_output"
        )

    reader = csv.DictReader(io.StringIO("\n".join(lines)), delimiter="\t")
    try:
        row = next(reader)
    except StopIteration:
        return SISTRResult(
            serovar=None, antigenic_formula=None, serogroup=None,
            h1=None, h2=None, o_antigen=None, error="parse_error"
        )

    return SISTRResult(
        serovar=clean_serovar(row.get("serovar")),
        antigenic_formula=row.get("antigenic_formula"),
        serogroup=row.get("serogroup"),
        h1=row.get("h1"),
        h2=row.get("h2"),
        o_antigen=row.get("o_antigen"),
        error=None,
    )
