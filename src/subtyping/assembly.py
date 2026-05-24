"""Download a single assembly FASTA from NCBI.

NCBI's assembly FTP layout:
    https://ftp.ncbi.nlm.nih.gov/genomes/all/{GCA|GCF}/{NNN}/{NNN}/{NNN}/{accession}_{name}/
        {accession}_{name}_genomic.fna.gz

We need to:
  1. Construct the directory URL from the accession (the prefix-splitting
     pattern is fixed)
  2. List the directory to find the actual filename (it includes the
     assembly name, not just the accession)
  3. Download the .fna.gz, decompress it, return the FASTA path

Assemblies are cached under cache/assemblies/{accession}.fna so we never
download the same one twice.
"""

from __future__ import annotations

import gzip
import logging
import re
import shutil
from pathlib import Path

import requests

from .. import config

log = logging.getLogger(__name__)


ASM_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all"


def _asm_dir_url(asm_acc: str) -> str:
    """Construct the NCBI assembly directory URL prefix from an accession.

    Example:
        GCA_000196035.1 → /genomes/all/GCA/000/196/035/
    """
    # GCA_000196035.1 → prefix GCA, then 000 196 035
    m = re.match(r"(GCA|GCF)_(\d{3})(\d{3})(\d{3})\.\d+", asm_acc)
    if not m:
        raise ValueError(f"Unrecognized assembly accession: {asm_acc}")
    prefix, a, b, c = m.groups()
    return f"{ASM_BASE}/{prefix}/{a}/{b}/{c}/"


def _find_fna_url(asm_acc: str) -> str | None:
    """List the NCBI assembly directory and find the .fna.gz file."""
    dir_url = _asm_dir_url(asm_acc)
    try:
        r = requests.get(dir_url, timeout=60)
        r.raise_for_status()
    except Exception as e:
        log.warning("Could not list %s: %s", dir_url, e)
        return None

    # Find the subdirectory matching the accession
    # The dir contains entries like "GCA_000196035.1_ASM19603v1/"
    m = re.search(rf'href="({re.escape(asm_acc)}[^"]*)/"', r.text)
    if not m:
        log.warning("No subdir found for %s in %s", asm_acc, dir_url)
        return None
    subdir = m.group(1)

    return f"{dir_url}{subdir}/{subdir}_genomic.fna.gz"


def fetch_assembly_fasta(asm_acc: str, cache_dir: Path | None = None) -> Path | None:
    """Download and decompress an assembly FASTA.

    Returns the path to the decompressed .fna file, or None on failure.
    Cached: if the file already exists locally, returns it without re-fetching.
    """
    cache_dir = cache_dir or (config.CACHE_DIR / "assemblies")
    cache_dir.mkdir(parents=True, exist_ok=True)
    fna_path = cache_dir / f"{asm_acc}.fna"
    if fna_path.exists() and fna_path.stat().st_size > 0:
        log.debug("Cached assembly: %s", fna_path)
        return fna_path

    gz_url = _find_fna_url(asm_acc)
    if not gz_url:
        return None

    gz_path = cache_dir / f"{asm_acc}.fna.gz"
    log.info("Downloading assembly: %s", gz_url)
    try:
        with requests.get(gz_url, stream=True, timeout=300) as r:
            if r.status_code == 404:
                log.warning("Assembly 404: %s", gz_url)
                return None
            r.raise_for_status()
            with open(gz_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=64 * 1024):
                    f.write(chunk)
    except Exception as e:
        log.warning("Assembly download failed for %s: %s", asm_acc, e)
        return None

    # Decompress
    try:
        with gzip.open(gz_path, "rb") as gz, open(fna_path, "wb") as out:
            shutil.copyfileobj(gz, out)
        gz_path.unlink()
    except Exception as e:
        log.warning("Decompression failed for %s: %s", asm_acc, e)
        if fna_path.exists():
            fna_path.unlink()
        return None

    return fna_path
