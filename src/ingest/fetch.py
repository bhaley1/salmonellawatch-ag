"""Download NCBI Pathogen Detection files.

Fetches three files per pathogen:
  - metadata.tsv        (canonical isolate metadata)
  - cluster_list.tsv    (PDT_acc → PDS_acc mapping)
  - amr.metadata.tsv    (AMR gene calls)

Caches downloads under cache/{taxgroup}/ with a freshness check.
NCBI updates files in place under the same release accession throughout
the day, so we re-download when the local copy is older than
CACHE_MAX_AGE_HOURS.
"""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass
from pathlib import Path

import requests

from .. import config

log = logging.getLogger(__name__)

PDG_RE = re.compile(r"PDG\d+\.\d+")


@dataclass
class Snapshot:
    """Files fetched for one pathogen on one run."""
    pathogen: str
    taxgroup: str
    pdg_release: str
    metadata_path: Path
    metadata_url: str
    metadata_bytes: int
    cluster_list_path: Path | None
    cluster_list_url: str | None
    cluster_list_bytes: int | None
    amr_path: Path | None
    amr_url: str | None
    amr_bytes: int | None


def _is_fresh(path: Path) -> bool:
    if not (path.exists() and path.stat().st_size > 0):
        return False
    age_seconds = time.time() - path.stat().st_mtime
    return age_seconds < (config.CACHE_MAX_AGE_HOURS * 3600)


def _latest_release(taxgroup: str) -> str | None:
    """Scrape NCBI's directory listing for the highest-numbered PDG release."""
    url = f"{config.NCBI_BASE}/{taxgroup}/"
    log.debug("Listing %s", url)
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
    except Exception as e:
        log.error("[%s] failed to list releases: %s", taxgroup, e)
        return None
    matches = PDG_RE.findall(r.text)
    if not matches:
        log.error("[%s] no PDG releases at %s", taxgroup, url)
        return None
    return sorted(
        set(matches),
        key=lambda s: tuple(int(x) for x in s[3:].split(".")),
    )[-1]


def _download(url: str, dest: Path, required: bool = True) -> Path | None:
    """Stream-download a file. If required=False, 404 returns None."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with requests.get(url, stream=True, timeout=600) as r:
            if r.status_code == 404:
                if required:
                    raise FileNotFoundError(url)
                log.warning("Optional file 404: %s", url)
                return None
            r.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 256):
                    f.write(chunk)
        return dest
    except Exception as e:
        if required:
            raise
        log.warning("Optional file fetch failed (%s): %s", url, e)
        return None


def fetch_one(pathogen: str) -> Snapshot | None:
    """Fetch all three files for one pathogen."""
    pcfg = config.PATHOGENS[pathogen]
    taxgroup = pcfg.taxgroup

    release = _latest_release(taxgroup)
    if not release:
        log.error("[%s] no release available; skipping", pathogen)
        return None
    log.info("[%s] latest release: %s", pathogen, release)

    tax_cache = config.CACHE_DIR / taxgroup
    snps_subdir = pcfg.snps_subdir if hasattr(pcfg, "snps_subdir") and pcfg.snps_subdir else None
    base_url = f"{config.NCBI_BASE}/{taxgroup}/{release}"
    snps_url = f"{config.NCBI_BASE}/{taxgroup}/{snps_subdir}" if snps_subdir else base_url

    # metadata.tsv
    md_path = tax_cache / f"{release}.metadata.tsv"
    md_url = f"{base_url}/Metadata/{release}.metadata.tsv"
    if _is_fresh(md_path):
        log.info("[%s] using cached metadata (fresh)", pathogen)
    else:
        if md_path.exists():
            log.info("[%s] metadata cache stale; re-fetching", pathogen)
        log.info("[%s] downloading %s", pathogen, md_url)
        _download(md_url, md_path, required=True)
    md_bytes = md_path.stat().st_size

    # cluster_list.tsv — try the canonical name first, fall back gracefully
    cl_path: Path | None = tax_cache / f"{release}.cluster_list.tsv"
    cl_url: str | None = None
    cl_bytes: int | None = None
    if _is_fresh(cl_path):
        log.info("[%s] using cached cluster_list (fresh)", pathogen)
        cl_url = f"{base_url}/Clusters/{release}.reference_target.cluster_list.tsv"
        cl_bytes = cl_path.stat().st_size
    else:
        # For pathogens with snps_subdir, cluster_list is under latest_snps
        # and uses a different release number - discover it from the directory
        if snps_subdir:
            import re as _re
            snps_dir_html = __import__("requests").get(f"{config.NCBI_BASE}/{taxgroup}/{snps_subdir}/Clusters/", timeout=30).text
            snps_releases = _re.findall(r"PDG\d+\.\d+", snps_dir_html)
            snps_rel = sorted(set(snps_releases), key=lambda s: tuple(int(x) for x in s[3:].split(".")))[-1] if snps_releases else release
        else:
            snps_rel = release
        candidates = [
            f"{snps_url}/Clusters/{snps_rel}.reference_target.cluster_list.tsv",
            f"{base_url}/Clusters/{release}.reference_target.cluster_list.tsv",
            f"{base_url}/Clusters/{release}.cluster_list.tsv",
        ]
        for url in candidates:
            got = _download(url, cl_path, required=False)
            if got and got.stat().st_size > 0:
                cl_url = url
                cl_bytes = got.stat().st_size
                log.info("[%s] cluster_list downloaded (%d bytes)", pathogen, cl_bytes)
                break
        else:
            log.warning("[%s] no cluster_list file available", pathogen)
            cl_path = None

    # amr.metadata.tsv
    amr_path: Path | None = tax_cache / f"{release}.amr.metadata.tsv"
    amr_url: str | None = None
    amr_bytes: int | None = None
    if _is_fresh(amr_path):
        log.info("[%s] using cached AMR (fresh)", pathogen)
        amr_url = f"{base_url}/AMR/{release}.amr.metadata.tsv"
        amr_bytes = amr_path.stat().st_size
    else:
        candidates = [
            f"{snps_url}/AMR/{snps_rel}.amr.metadata.tsv",
            f"{base_url}/AMR/{release}.amr.metadata.tsv",
            f"{base_url}/AMRFinderPlus/{release}.amr.metadata.tsv",
        ]
        for url in candidates:
            got = _download(url, amr_path, required=False)
            if got and got.stat().st_size > 0:
                amr_url = url
                amr_bytes = got.stat().st_size
                log.info("[%s] AMR downloaded (%d bytes)", pathogen, amr_bytes)
                break
        else:
            log.warning("[%s] no AMR file available", pathogen)
            amr_path = None

    return Snapshot(
        pathogen=pathogen,
        taxgroup=taxgroup,
        pdg_release=release,
        metadata_path=md_path,
        metadata_url=md_url,
        metadata_bytes=md_bytes,
        cluster_list_path=cl_path,
        cluster_list_url=cl_url,
        cluster_list_bytes=cl_bytes,
        amr_path=amr_path,
        amr_url=amr_url,
        amr_bytes=amr_bytes,
    )


def fetch_all() -> list[Snapshot]:
    out: list[Snapshot] = []
    for pathogen in config.PATHOGENS:
        snap = fetch_one(pathogen)
        if snap:
            out.append(snap)
        time.sleep(1)  # be polite to NCBI
    return out
