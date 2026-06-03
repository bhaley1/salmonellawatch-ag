# SalmonellaWatch-Ag

A daily-updated digest of *Salmonella enterica* genomic clusters from
[NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/) with
at least one non-human (food, animal, or environmental) isolate collected
in the last 180 days. Designed for agricultural food safety and
pre-harvest surveillance.

**Live dashboard:** https://bhaley1.github.io/salmonellawatch-ag/

**How-to guide:** https://bhaley1.github.io/salmonellawatch-ag/guide.html

**Methods & definitions:** https://bhaley1.github.io/salmonellawatch-ag/methods.html

## What this is

SalmonellaWatch-Ag is a personal research project that ingests the daily
NCBI Pathogen Detection metadata for *Salmonella enterica*, identifies
SNP clusters with recent non-human isolates, and surfaces them on a
dashboard alongside signals that may indicate clusters worth investigating:
geographic spread, possible imported food vehicles, travel-associated
cases, multi-submitter detection, AMR genotypes, acceleration in deposit
rate, emerging human cases in a previously food-dominant cluster, MLST
sequence typing, and commodity-coded source profiles.

The agricultural focus means the dashboard surfaces clusters relevant to
pre-harvest food safety, environmental monitoring, and commodity-source
tracking, including clusters that may not yet have associated human illness.

Key features:

- ~515 active clusters with at least one non-human isolate in the last 180 days
- Filterable by country, serotype, MLST sequence type, and Serotypes of Concern (Katz et al. 2024)
- Commodity-coded source profiles and stacked-bar histograms (poultry, bovine, swine, produce, surface water, human, other)
- Signal chips flagging clusters that are accelerating, internationally distributed, transitioning from food/environment to human cases, or associated with imported food or travel
- Geographic footprint maps and per-isolate detail with links to NCBI BioSample records
- MLST typing (Achtman 7-gene scheme) for 497 of 500 active clusters

## What this is NOT

- **Not an official surveillance system.** USDA-FSIS, CDC, and FDA issue
  health advisories, food recalls, and outbreak declarations. This
  dashboard does not.
- **Not a regulatory decision-support tool.** Genomic cluster membership
  does not prove a foodborne link. Confirming an outbreak requires
  epidemiological investigation.
- **Not a peer-reviewed product.** The signal thresholds, commodity
  classifications, and SOC designations reflect the author's compilation
  of public literature, not peer-reviewed claims.

See the [methods page](https://bhaley1.github.io/salmonellawatch-ag/methods.html)
for full documentation of every signal, definition, threshold, and caveat.

## Relationship to SalmonellaWatch

SalmonellaWatch-Ag is an agricultural mirror of [SalmonellaWatch](https://bhaley1.github.io/salmonellawatch/),
which focuses on clusters with recent human cases (60-day window). Both
sites share the same underlying SQLite database and rendering pipeline.

| | SalmonellaWatch | SalmonellaWatch-Ag |
|---|---|---|
| Focus | Human public health | Agricultural food safety |
| Active cluster definition | >= 1 human case in last 60 days | >= 1 non-human isolate in last 180 days |
| Primary audience | Epidemiologists, public health labs | Food safety microbiologists, pre-harvest researchers |
| Commodity histograms | No | Yes |
| Source profile badges | No | Yes |

## Data sources

All cluster, isolate, and AMR data are sourced from
[NCBI Pathogen Detection](https://www.ncbi.nlm.nih.gov/pathogens/).
Three files per release:

- `metadata.tsv` -- per-isolate metadata (host, source, geography, dates, IFSAC category)
- `cluster_list.tsv` -- per-cluster isolate membership
- `amr.metadata.tsv` -- AMRFinderPlus gene calls per isolate

## Architecture

Single SQLite database. Three pipelines: ingest, subtyping, render.
NCBI FTP
|
ingest         metadata.tsv -> isolates table
cluster_list.tsv -> cluster membership
amr.metadata.tsv -> AMR gene calls
|
subtyping      mlst (Torsten Seemann, Achtman scheme) -> ST
representative isolate picked per cluster
cluster_typing table
|
summarize      per-cluster signals computed
commodity classification (keyword-based)
cluster_summary table
|
render         Jinja2 -> static HTML
site/ (deployed to GitHub Pages via docs/)

## Repository structure
src/
config.py              # tuning knobs (RECENT_WINDOW_DAYS = 180)
signals.py             # signal functions (acceleration, emergence, etc.)
lookups/               # geography, AMR, centroids
ingest/                # NCBI fetch/parse/upsert/summarize
subtyping/             # MLST typing
render/
queries.py           # SQLite queries, commodity classification, histogram builder
run.py               # Jinja2 rendering
map.py               # per-cluster SVG world maps
templates/
index.html             # main dashboard
methods.html           # methods & definitions
guide.html             # how-to guide for microbiologists
assets/css/main.css    # all styling
docs/                    # GitHub Pages deployment (copy of site/)
cache/
assemblies/            # downloaded NCBI assemblies for MLST
mlst_results.json      # cached MLST calls
db/                      # SQLite database (symlinked, not in repo)

## Running locally

```bash
# Activate conda environment
conda activate listeriawatch

# Full pipeline
python -m src.ingest.run -v
python -m src.subtyping.run -v
python -m src.render.run -v

# Preview
cd site && python3 -m http.server 8002
# open http://localhost:8002
```

## Commodity classification

Non-human isolates are classified into commodity categories by keyword
matching against `isolation_source` and `ifsac_category` fields:

| Category | Keywords |
|---|---|
| Poultry | chicken, turkey, poultry, broiler, avian, duck, goose, hen, egg |
| Bovine | beef, bovine, cattle, cow, calf, veal, dairy, ruminant |
| Swine | pork, swine, pig, sow, boar, piglet, hog |
| Produce | vegetable, fruit, leafy, lettuce, spinach, tomato, sprout, herb, nut, melon |
| Water | water, aquatic, stream, river, lake, pond |
| Other | all other non-human sources |

## Serotypes of Concern (SOC)

The 21 SOC serotypes follow the framework of Katz et al. (2024):

> Katz TS, Harhay DM, Schmidt JW, Wheeler TL. Identifying a list of
> Salmonella serotypes of concern to target for reducing risk of
> salmonellosis. *Front Microbiol.* 2024 Feb 12;15:1307563.
> [doi:10.3389/fmicb.2024.1307563](https://doi.org/10.3389/fmicb.2024.1307563)

SOC serotypes: Enteritidis, Typhimurium, I 4,[5],12:i:- (monophasic),
Heidelberg, Infantis, Newport, Uganda, Braenderup, Muenchen, Montevideo,
Javiana, Reading, Dublin, Oranienburg, Potsdam, Thompson, Saintpaul,
Hadar, Schwarzengrund, Anatum, Berta.

## License

The code in this repository is MIT-licensed. The dashboard outputs
include public-domain NCBI data; downstream display of that data is
subject to NCBI's terms of use.

## Acknowledgments

- NCBI Pathogen Detection team for the data
- Torsten Seemann's [mlst](https://github.com/tseemann/mlst) tool
- AMRFinderPlus (NCBI)
- PubMLST for the *Salmonella enterica* MLST scheme (Achtman)
- Katz et al. (2024) for the SOC framework
- All the public-health and food-safety labs whose deposits make this possible
