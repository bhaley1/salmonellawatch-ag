-- Pathogen Watch v2 schema
-- One-shot canonical schema. Run with sqlite3 to initialize.
-- All tables use SQLite's INTEGER PRIMARY KEY ROWID semantics where appropriate.

-- ---------------------------------------------------------------------------
-- isolates: one row per PDT_acc per pathogen
-- ---------------------------------------------------------------------------
-- target_creation_date is the column we use for "recent activity":
-- it's NCBI's timestamp for when the isolate first entered Pathogen Detection,
-- which is the right semantic for "this case became visible to surveillance on date X".
-- collection_date is when the sample was actually taken, which can be months earlier.
CREATE TABLE IF NOT EXISTS isolates (
    pdt_acc                 TEXT PRIMARY KEY,
    pathogen                TEXT NOT NULL,
    pds_acc                 TEXT,
    epi_type                TEXT,              -- 'clinical' | 'environmental/other' | NULL
    source_category         TEXT,              -- derived: 'Human' | 'Food' | 'Animal' | 'Environment' | 'Unknown'
    host                    TEXT,
    isolation_source        TEXT,
    geo_loc_name            TEXT,              -- raw NCBI string
    geo_country             TEXT,              -- parsed: 'USA', 'United Kingdom', etc.
    geo_admin1              TEXT,              -- parsed: 'Maryland', 'Bavaria', etc., NULL if not present
    collection_date         DATE,              -- parsed; NULL if unparseable
    collection_date_raw     TEXT,              -- preserve original NCBI string
    target_creation_date    DATE,              -- when NCBI ingested this isolate
    scientific_name         TEXT,
    serovar                 TEXT,
    biosample_acc           TEXT,
    asm_acc                 TEXT,
    sra_acc                 TEXT,
    asm_level               TEXT,              -- 'Complete Genome' | 'Chromosome' | 'Scaffold' | 'Contig' | NULL
    asm_stats_contig_n50    INTEGER,
    food_origin             TEXT,              -- food's country of origin (when populated by submitter)
    ifsac_category          TEXT,              -- IFSAC food category
    host_disease            TEXT,              -- NCBI host_disease field; sometimes contains travel notes
    bioproject_acc          TEXT,              -- BioProject accession (PRJNA…); submitter diversity proxy
    pdg_release             TEXT NOT NULL,     -- which NCBI release we last saw this in
    last_seen_at            TIMESTAMP NOT NULL -- updated every ingest run
);

CREATE INDEX IF NOT EXISTS idx_isolates_pds_acc ON isolates(pds_acc);
CREATE INDEX IF NOT EXISTS idx_isolates_pathogen ON isolates(pathogen);
CREATE INDEX IF NOT EXISTS idx_isolates_target_creation_date ON isolates(target_creation_date);
CREATE INDEX IF NOT EXISTS idx_isolates_epi_type ON isolates(epi_type);
CREATE INDEX IF NOT EXISTS idx_isolates_pathogen_epi ON isolates(pathogen, epi_type);


-- ---------------------------------------------------------------------------
-- isolate_amr: per-isolate AMR gene calls from NCBI's AMR sidecar file
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS isolate_amr (
    pdt_acc                 TEXT NOT NULL,
    gene_symbol             TEXT NOT NULL,
    element_subtype         TEXT,              -- e.g., 'AMR', 'STRESS', 'VIRULENCE'
    PRIMARY KEY (pdt_acc, gene_symbol)
);

CREATE INDEX IF NOT EXISTS idx_isolate_amr_pdt ON isolate_amr(pdt_acc);


-- ---------------------------------------------------------------------------
-- cluster_typing: cluster-level consensus typing
-- ---------------------------------------------------------------------------
-- Populated by the typing pipeline (Phase 2).
-- Persists across releases because (pathogen, pds_acc) is stable.
CREATE TABLE IF NOT EXISTS cluster_typing (
    pathogen                TEXT NOT NULL,
    pds_acc                 TEXT NOT NULL,

    consensus_serovar       TEXT,
    consensus_serovar_n     INTEGER,           -- how many cluster members agreed
    consensus_serovar_total INTEGER,           -- how many had any non-null serovar

    mlst_scheme             TEXT,
    mlst_st                 TEXT,              -- the ST call: 'ST6', 'novel', 'untypeable', etc.
    mlst_alleles            TEXT,              -- JSON: {gene: allele_number}
    mlst_representative_pdt TEXT,              -- which isolate we typed
    mlst_error              TEXT,              -- 'no_assembly' | 'no_match' | NULL

    typed_at                TIMESTAMP NOT NULL,
    PRIMARY KEY (pathogen, pds_acc)
);


-- ---------------------------------------------------------------------------
-- cluster_summary: materialized per-cluster summary
-- ---------------------------------------------------------------------------
-- Rebuilt after every ingest. Lets the renderer query in one shot.
CREATE TABLE IF NOT EXISTS cluster_summary (
    pathogen                       TEXT NOT NULL,
    pds_acc                        TEXT NOT NULL,

    n_total                        INTEGER NOT NULL,
    n_human                        INTEGER NOT NULL,
    n_nonhuman                     INTEGER NOT NULL,
    n_food                         INTEGER NOT NULL DEFAULT 0,
    n_animal                       INTEGER NOT NULL DEFAULT 0,
    n_environment                  INTEGER NOT NULL DEFAULT 0,

    earliest_collection_date       DATE,
    latest_collection_date         DATE,
    earliest_target_creation_date  DATE,
    latest_target_creation_date    DATE,
    temporal_span_days             INTEGER,   -- latest_collection_date - earliest_collection_date, NULL if missing dates

    -- Oldest-isolate signatures: the "when did this cluster first appear" details.
    -- We store three separately because the epidemiological meaning differs:
    --   oldest of any kind = when did the signal first appear at all
    --   oldest human       = when did the cluster first cause disease
    --   oldest nonhuman    = when did the (potential) source first show up
    oldest_isolate_json            TEXT,      -- {pdt, date, geo, source, source_category}
    oldest_human_json              TEXT,      -- {pdt, date, geo, source}
    oldest_nonhuman_json           TEXT,      -- {pdt, date, geo, source, source_category}

    -- Histogram: per-year counts of human vs. nonhuman isolates
    histogram_json                 TEXT,      -- [{year, n_human, n_nonhuman}, ...]
    histogram_max_year_count       INTEGER,   -- biggest single-year total, for chart scaling

    -- Collection-to-deposit lag (in days) summarizing how fresh the cluster's
    -- data is. Median is the headline; mean is also stored for completeness.
    deposit_lag_median             INTEGER,
    deposit_lag_mean               INTEGER,

    -- Animal host species breakdown — surfaces actual host strings (Bos taurus,
    -- Sus scrofa, etc.) that would otherwise be hidden behind the 'Animal'
    -- source category aggregation. Only populated when nonhuman hosts exist.
    host_summary_json              TEXT,

    -- Per-cluster geographic footprint over the last 365 days.
    -- Aggregated by US state (when admin1 available) or country.
    -- Pre-rendered as a list of split-dot location dicts.
    map_locations_json             TEXT,

    -- Most-recently-deposited isolate in this cluster that has an assembly
    -- accession. Used to show users a "view the latest assembled genome
    -- for this cluster" link.
    latest_assembly_json           TEXT,      -- {pdt, biosample, asm_acc, collection_date, target_creation_date}

    -- Per-cluster derived signals: geographic spread, import signal, AMR
    -- critical resistance, acceleration, IFSAC summary, clonal complex.
    -- All computed in signals.py and stored as a single JSON blob.
    signals_json                   TEXT,

    countries_json                 TEXT,       -- JSON array of country strings, sorted by count desc
    admin1_json                    TEXT,       -- JSON: {by_country: {country: [{admin1, n}]}, unspecified: {country: n}}
    source_summary_json            TEXT,       -- JSON array of {category, source, n}

    new_humans_in_window           INTEGER NOT NULL DEFAULT 0,
    new_humans_30d                 INTEGER NOT NULL DEFAULT 0,
    new_humans_15d                 INTEGER NOT NULL DEFAULT 0,
    new_humans_in_window_pdts_json  TEXT,       -- JSON array of pdt_accs
    new_humans_in_window_dates_json TEXT,       -- JSON array of {pdt, date_added, collection, geo}
    window_days                    INTEGER NOT NULL DEFAULT 60,  -- the window this row reflects

    refreshed_at                   TIMESTAMP NOT NULL,
    PRIMARY KEY (pathogen, pds_acc)
);

CREATE INDEX IF NOT EXISTS idx_cluster_summary_new_humans ON cluster_summary(new_humans_in_window DESC);
CREATE INDEX IF NOT EXISTS idx_cluster_summary_pathogen ON cluster_summary(pathogen);


-- ---------------------------------------------------------------------------
-- pdg_releases: which NCBI releases we've processed and when
-- ---------------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS pdg_releases (
    pathogen                TEXT NOT NULL,
    pdg_release             TEXT NOT NULL,
    metadata_url            TEXT,
    metadata_bytes          INTEGER,
    cluster_list_url        TEXT,
    cluster_list_bytes      INTEGER,
    amr_url                 TEXT,
    amr_bytes               INTEGER,
    ingested_at             TIMESTAMP NOT NULL,
    PRIMARY KEY (pathogen, pdg_release)
);
