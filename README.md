# ListeriaWatch

A daily-updated digest of *Listeria monocytogenes* genomic clusters from
[NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/) that
have at least one human case collected in the last 60 days.

**Live dashboard:** https://bhaley1.github.io/listeria-watch/
**Methods & definitions:** https://bhaley1.github.io/listeria-watch/methods.html

## What this is

ListeriaWatch is a personal research project that ingests the daily
NCBI Pathogen Detection metadata for *Listeria monocytogenes*,
identifies SNP clusters with recent human cases, and surfaces them on
a single dashboard alongside signals that may indicate clusters worth
investigating: geographic spread, possible imported food vehicles,
travel-associated cases, multi-submitter detection, AMR gene presence,
acceleration in deposit rate, emerging human cases in a previously
food-dominant cluster, and MLST/CC typing.

The dashboard updates nightly via a GitHub Actions cron job.

## What this is NOT

- **Not an official surveillance system.** CDC, FDA, and USDA-FSIS issue
  health advisories, food recalls, and outbreak declarations. This
  dashboard does not.
- **Not a clinical decision-support tool.** AMR gene presence on the
  dashboard does not equal confirmed phenotypic resistance.
- **Not a peer-reviewed product.** The CC labels, AMR gene curation,
  and signal thresholds reflect the author's compilation of public
  literature, not peer-reviewed claims.

See the [methods page](https://bhaley1.github.io/listeria-watch/methods.html)
for full documentation of every signal, definition, threshold, and
caveat.

## Architecture

Single SQLite database. Three pipelines: ingest, subtyping, render.

```
NCBI FTP
   ↓
ingest         metadata.tsv → isolates table
              cluster_list.tsv → cluster membership
              amr.metadata.tsv → AMR gene calls
   ↓
subtyping     mlst (Torsten Seemann) → ST → CC lookup
              representative isolate picked per cluster
              cluster_typing table
   ↓
summarize     per-cluster signals computed
              cluster_summary table
   ↓
render        Jinja2 → static HTML
              site/ (deployed to GitHub Pages)
```

Tests: 127 across `tests/test_pipeline.py`, `tests/test_signals.py`,
`tests/test_subtyping.py`, `tests/test_map.py`. All synthetic-fixture
based; real NCBI integration is exercised manually.

## Repository structure

```
src/
  config.py              # all tuning knobs
  schema.sql             # SQLite schema
  signals.py             # the nine signal functions
  lookups/               # geography, AMR, CC, centroids
  ingest/                # NCBI fetch/parse/upsert/summarize
  subtyping/             # MLST and serovar typing
  render/                # Jinja templates → static site
templates/
  index.html             # main dashboard
  methods.html           # methods & definitions
  assets/css/main.css    # all styling
tests/
  test_*.py              # 127 tests
  fixtures/              # synthetic Listeria fixtures
.github/workflows/
  daily.yml              # nightly cron build + deploy
```

## Running locally

```bash
# Setup
python -m pip install -r requirements.txt
# Need mlst (and BLAST+) for full typing:
# https://github.com/tseemann/mlst

# Full pipeline
python -m src.ingest.run -v
python -m src.subtyping.run -v
python -m src.render.run -v

# Preview
cd site && python3 -m http.server 8000
# open http://localhost:8000
```

## Review mode

When `LW_REVIEW_MODE=1` (the default), the dashboard shows a "PRIVATE
REVIEW — DO NOT SHARE" banner at the top of every page. Set it to `0`
when ready to go public:

- Locally: `LW_REVIEW_MODE=0 python -m src.render.run`
- In the workflow: edit `.github/workflows/daily.yml`

## License

The code in this repository is MIT-licensed. The dashboard outputs
include public-domain NCBI data; downstream display of that data is
subject to NCBI's terms of use.

## Acknowledgments

- NCBI Pathogen Detection team for the data
- Torsten Seemann's [mlst](https://github.com/tseemann/mlst) tool
- AMRFinderPlus (NCBI)
- PubMLST for the *Listeria monocytogenes* MLST scheme
- All the public-health labs whose deposits make this possible
