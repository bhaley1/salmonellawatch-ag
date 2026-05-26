# Deployment checklist — ListeriaWatch v2

Walk through this once to deploy v2 to GitHub Pages alongside the
existing v1 (which stays running until you validate v2).

Estimated time: **2–3 hours** including the local smoke test.

---

## Phase A — Local smoke test against real NCBI data

You haven't run phase 8 against live NCBI since several rounds of
changes. Do this FIRST, before any deployment. If something breaks on
real data, we want to find out privately.

### A1. Pull the latest v2 code

```bash
cd /Users/braddhaley/Downloads/pw-v2/
# Unzip the latest deliverable if needed
```

### A2. Set up the environment

```bash
# Use the M3 Max miniforge3 setup
conda activate <your-listeriawatch-env>   # or create one
pip install -r requirements.txt

# Verify mlst is installed (you set this up in earlier sessions)
which mlst
mlst --version
```

### A3. Run the full pipeline

```bash
# Clear any stale DB from earlier rounds
rm -f db/pathogen-watch.sqlite

# Real ingest from NCBI
python -m src.ingest.run -v
# Expected: ~20-40 minutes; downloads metadata.tsv, cluster_list.tsv,
# amr.metadata.tsv for Listeria. Logs should show "Found N clusters"
# and "Upserted M isolates".

# Typing (MLST + CC)
python -m src.subtyping.run -v --mlst-budget 500
# Expected: ~30-60 minutes if many new clusters need typing.
# Use Ctrl+C if it's taking too long — the budget cap prevents runaway.

# Summarization (fast)
python -c "
from src import db
from src.ingest import summarize
import datetime
with db.connect() as conn:
    summarize.materialize_cluster_summary(conn, today=datetime.date.today())
"

# Render
python -m src.render.run -v
```

### A4. Preview locally

```bash
cd site
python3 -m http.server 8000
# Open http://localhost:8000 in your browser
# Open http://localhost:8000/methods.html
```

### A5. What to look for

- [ ] Active cluster count is reasonable (probably 10-50 for Listeria)
- [ ] Each cluster card renders without obviously broken sections
- [ ] The world map shows continent silhouettes with dots placed correctly
- [ ] Country dropdown filter populates with real countries
- [ ] Signal chips appear in reasonable patterns:
  - MULTI-STATE / MULTI-COUNTRY on clusters that genuinely span borders
  - IMPORT SIGNAL only when food_origin or country mismatch is real
  - ACCELERATING only on clusters with rapid recent deposits
  - EMERGING only on clusters with the food→human transition pattern
  - MULTI-SUBMITTER only on clusters with ≥5 BioProjects
  - CC labels only on CC1/CC2/CC4/CC6/CC9/CC121
- [ ] PDT accession links open NCBI BioSample pages in new tabs
- [ ] "Latest assembly" link goes to a valid NCBI Datasets page
- [ ] Methods page renders fully with all sections
- [ ] References section links all work

### A6. If something breaks

Send a screenshot. We'll iterate before any public deploy.

---

## Phase B — Create the new GitHub repository

### B1. Create the repo

On github.com:
- Repository name: `listeria-watch`
- Owner: `bhaley1`
- Public visibility
- Do not initialize with a README (we have one)
- Do not add .gitignore or license (we'll add them ourselves if needed)

### B2. Push the code

```bash
cd /Users/braddhaley/Downloads/pw-v2/
git init
git add .
git commit -m "Initial commit: ListeriaWatch v2"
git branch -M main
git remote add origin git@github.com:bhaley1/listeria-watch.git
git push -u origin main
```

### B3. Enable GitHub Pages

On the repo's Settings → Pages:
- Source: Deploy from a branch
- Branch: `gh-pages` (will not exist yet — that's fine)
- Folder: `/ (root)`
- Click Save

After the first successful Actions run (Phase C), the branch will exist
and the site will go live at `https://bhaley1.github.io/listeria-watch/`.

---

## Phase C — Trigger the first deploy

### C1. Manually trigger the workflow

On the repo's Actions tab:
- Find "Daily build" workflow
- Click "Run workflow" → "Run workflow"

### C2. Watch the run

Expect:
- "Set up Python" — ~30s
- "Install Python dependencies" — ~30s
- "Install BLAST+" — ~30s
- "Install mlst" — ~10s
- "Restore caches" — ~5s (cold cache, fast)
- "Run ingest" — **~20-40 minutes** (downloading from NCBI)
- "Run subtyping" — **~10-20 minutes** (with --mlst-budget 200)
- "Run summarization" — ~30s
- "Render site" — ~30s
- "Deploy to gh-pages" — ~30s

Total: ~30-60 minutes for the first run, ~5-15 minutes for subsequent
cached runs.

### C3. Verify the site is live

Wait 1-2 minutes after deploy step completes, then visit:
- https://bhaley1.github.io/listeria-watch/
- https://bhaley1.github.io/listeria-watch/methods.html

You should see:
- The yellow "PRIVATE REVIEW — DO NOT SHARE" banner at the top
- Below that the standard "Beta" banner
- The full dashboard with real Listeria clusters

---

## Phase D — Share with colleagues

The URL is `https://bhaley1.github.io/listeria-watch/`.

### D1. Pick reviewers

Send the URL only to people who'll give feedback. The review banner
tells them not to circulate it. A small number (3-8 people) is
probably right for the first review round.

### D2. Suggested feedback ask

Something like:

> I've built a daily-updated dashboard for *Listeria monocytogenes*
> SNP clusters from NCBI Pathogen Detection. It surfaces clusters
> with recent human cases and flags signals worth investigating:
> geographic spread, possible imports, acceleration, AMR, emerging
> human cases in food-dominant clusters, etc.
>
> Before public release I'd like your eye on:
>
> 1. **CC labels** (CC1, CC2, CC4, CC6, CC9, CC121 — see methods page).
>    Are my one-line summaries accurate?
> 2. **AMR gene curation** (full list on methods page). Are these the
>    right first-line-treatment genes to flag for *Listeria*?
> 3. **Signal thresholds** — acceleration at 2×, emergence at 30%
>    human + 20pp jump, multi-submitter at 5+ BioProjects. Reasonable?
> 4. **Anything else that looks wrong**, including data interpretation,
>    epidemiology framing, methods-page caveats.
>
> URL: https://bhaley1.github.io/listeria-watch/
>
> Don't circulate the URL publicly — there's a banner explaining this
> is a private review.

---

## Phase E — Iterate based on feedback

Make changes locally. Commit and push to `main`. The workflow
auto-deploys on push.

### E1. Common iteration cycle

```bash
# Make changes
git add -A && git commit -m "fix: CC4 label per feedback from Smith"
git push
# Workflow runs in 2-5 minutes
```

### E2. When ready to go public

1. Edit `.github/workflows/daily.yml`
2. Change `LW_REVIEW_MODE: "1"` to `LW_REVIEW_MODE: "0"`
3. Commit and push
4. The next build removes the review banner

---

## Phase F — Relationship with v1

While v2 is in review, v1 at `https://bhaley1.github.io/pathogen-watch/`
keeps running. Do not touch v1 yet.

### F1. After v2 is validated and public

Three options for v1:

**Option F1a — Keep v1 indefinitely** at its current URL. Mention v2 as
a successor in v1's README. Some users may prefer v1's older signals.

**Option F1b — Redirect v1 to v2.** Add a small index.html in v1 that
redirects to v2. v1's URL keeps working.

**Option F1c — Retire v1.** Disable GitHub Pages on the v1 repo. URL
breaks (users get a 404). Only do this after several weeks of v2 being
stable.

My recommendation: **F1b (redirect)**. Lowest-friction for existing
users.

---

## Troubleshooting

### "The first cron run took 90 minutes and timed out"

The default GitHub Actions timeout is 6 hours, so this shouldn't
happen. If it does, lower `--mlst-budget` in `daily.yml` from 200 to
100 and re-run; the cached state will pick up where the previous run
left off.

### "The map is empty"

Check whether the cluster has any isolates with `geo_country` and
`collection_date` set. The map filter is "all isolates ever with a
location" — if everyone is in clusters whose members lack country or
collection date, the map will be empty.

### "MLST shows 'pending' for everything"

Check that mlst is in PATH on the Actions runner. The
`Install mlst` step prepends `/opt/mlst/bin` to `$GITHUB_PATH`; if
something interferes with that, mlst won't run.

### "I want to undo a bad commit on gh-pages"

The `force_orphan: true` flag in the deploy step means we don't keep
history on `gh-pages` — every deploy rewrites it. Just push a fix to
main; the next deploy will replace the bad content within 2-5 minutes.

### "v2 looks really different from v1 and a colleague asks why"

Point them to the methods page. v2's signal set is broader, the map
is real, and the cluster cards have more information per card. v1
predates the v2 work in this conversation.
