"""
SalmonellaWatch-Ag cluster matcher — Streamlit GUI

Drop one or more FASTA/FNA assemblies; the app runs the full typing panel on each
(cluster match, SISTR serovar, MLST, plasmid replicons, pESI screen, AMRFinderPlus)
and renders the result as a card matching SalmonellaWatch-Ag's visual language.

Run:
    streamlit run app.py

Configuration is read from settings.json next to this file. If absent, sensible
defaults are derived from the project layout.
"""

import json
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import streamlit as st

# ---------- paths and settings ----------

HERE = Path(__file__).parent
RESULTS_DIR = HERE / "results"
RESULTS_DIR.mkdir(exist_ok=True)
SETTINGS_PATH = HERE / "settings.json"

DEFAULT_SETTINGS = {
    "type_genome_script": str(HERE / "type_genome.py"),
    "snp_cluster_script": str(HERE / "snp_cluster_match.py"),
    "index_npz":        str(HERE / "cluster_index.npz"),
    "schema_frozen":    str(HERE / "cache" / "schema_frozen"),
    "core_loci":        str(HERE / "cache" / "cgmlst_extract" / "cgMLSTschema95.txt"),
    "pathogen_watch_db": "/Users/braddhaley/Downloads/salmonellawatch/db/pathogen-watch.sqlite",
    "salmonellawatch_url": "https://bhaley1.github.io/salmonellawatch-ag/",
    "ncbi_pds_url_pattern": "https://www.ncbi.nlm.nih.gov/pathogens/isolates/#{pds}",
    "index_build_date": "2026-06-10",   # update when the index is rebuilt
    "cpu": 14,
}

def load_settings():
    if SETTINGS_PATH.exists():
        try:
            with open(SETTINGS_PATH) as f:
                return {**DEFAULT_SETTINGS, **json.load(f)}
        except Exception:
            pass
    return DEFAULT_SETTINGS.copy()

S = load_settings()

# ---------- streamlit page setup ----------

st.set_page_config(
    page_title="SalmonellaWatch-Ag cluster matcher",
    page_icon="🧬",
    layout="wide",
)

# tiny CSS to match the dashboard's visual cues
st.markdown("""
<style>
.badge { display: inline-block; padding: 2px 10px; border-radius: 4px;
         font-size: 11px; font-weight: 600; letter-spacing: 0.4px;
         margin-right: 6px; margin-bottom: 4px; }
.badge-soc    { background: #fee2e2; color: #991b1b; border: 1px solid #fca5a5; }
.badge-signal { background: #fef3c7; color: #92400e; border: 1px solid #fcd34d; }
.badge-multi  { background: #dbeafe; color: #1e3a8a; border: 1px solid #93c5fd; }
.badge-accel  { background: #fce7f3; color: #9d174d; border: 1px solid #f9a8d4; }
.badge-ok     { background: #d1fae5; color: #065f46; border: 1px solid #6ee7b7; }
.method-EXACTX { color: #065f46; font-weight: 600; }
.method-BLASTX { color: #1e40af; }
.method-PARTIALX { color: #92400e; font-weight: 600; }
.method-INTERNAL_STOP, .method-POINTX { color: #991b1b; font-weight: 600; }
.cluster-id { font-family: ui-monospace, SFMono-Regular, monospace;
              background: #f3f4f6; padding: 1px 6px; border-radius: 3px; }
.caveat-box { background: #fef3c7; border-left: 4px solid #f59e0b;
              padding: 10px 14px; margin: 8px 0; border-radius: 3px; }
</style>
""", unsafe_allow_html=True)

# ---------- helpers ----------

def signal_badges(signals):
    icons = {
        "geographic_spread": ("MULTI HEMISPHERE", "badge-multi"),
        "import_signal":     ("IMPORT SIGNAL",    "badge-signal"),
        "submitter_diversity": ("MULTI-SUBMITTER", "badge-multi"),
        "acceleration":      ("ACCELERATING",     "badge-accel"),
        "amr_critical":      ("AMR CRITICAL",     "badge-soc"),
        "human_emergence":   ("HUMAN EMERGENCE",  "badge-soc"),
        "travel_signal":     ("TRAVEL SIGNAL",    "badge-signal"),
    }
    out = []
    for s in signals or []:
        if s in icons:
            label, cls = icons[s]
            out.append(f'<span class="badge {cls}">{label}</span>')
    return "".join(out)

def run_panel(assembly_path, cpu, status_cb=None):
    """Run type_genome.py on one assembly; return parsed dict or {'error': ...}."""
    cmd = [
        sys.executable, S["type_genome_script"],
        "--assembly", str(assembly_path),
        "--index", S["index_npz"],
        "--schema", S["schema_frozen"],
        "--core-loci", S["core_loci"],
        "--db", S["pathogen_watch_db"],
        "--cpu", str(cpu),
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            return {"error": f"type_genome.py failed: {proc.stderr[:500]}"}
        return json.loads(proc.stdout)
    except subprocess.TimeoutExpired:
        return {"error": "timed out (>10 min)"}
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"}

def save_result(name, record):
    day = datetime.now().strftime("%Y-%m-%d")
    day_dir = RESULTS_DIR / day
    day_dir.mkdir(exist_ok=True)
    stem = Path(name).stem
    ts = datetime.now().strftime("%H%M%S")
    out = day_dir / f"{stem}_{ts}.json"
    with open(out, "w") as f:
        json.dump(record, f, indent=2)
    return out

def list_history():
    rows = []
    for day_dir in sorted(RESULTS_DIR.iterdir(), reverse=True):
        if not day_dir.is_dir():
            continue
        for p in sorted(day_dir.iterdir(), reverse=True):
            if p.suffix == ".json":
                rows.append(p)
    return rows[:30]

# ---------- result rendering ----------

def render_cluster_section(c, settings):
    if c.get("status") != "ok":
        st.error(f"Cluster call failed: {c.get('error', c.get('reason', 'unknown'))}")
        return

    pds = c.get("pds", "(none)")
    enrich = c.get("enrich") or {}
    serovar = enrich.get("serovar", "?")
    mlst_from_cluster = enrich.get("mlst_st_cluster", "?")

    # confidence flags
    flag_bits = []
    if c.get("ambiguous"):       flag_bits.append('<span class="badge badge-signal">AMBIGUOUS</span>')
    if c.get("no_close_match"):  flag_bits.append('<span class="badge badge-soc">NO CLOSE MATCH</span>')
    if c.get("low_overlap"):     flag_bits.append('<span class="badge badge-soc">LOW OVERLAP</span>')
    if not flag_bits:            flag_bits.append('<span class="badge badge-ok">CONFIDENT</span>')

    st.markdown(
        f"### {serovar} · {mlst_from_cluster} · "
        f'<span class="cluster-id">{pds}</span>',
        unsafe_allow_html=True,
    )
    st.markdown(" ".join(flag_bits), unsafe_allow_html=True)

    col1, col2, col3 = st.columns(3)
    col1.metric("Allele distance", c.get("distance", "—"))
    col2.metric("Margin to next cluster", c.get("margin", "—"))
    col3.metric("Next cluster", c.get("next_pds", "—"))

    # deep links
    ncbi_url = settings["ncbi_pds_url_pattern"].format(pds=pds)
    dashboard_url = settings["salmonellawatch_url"]
    lc1, lc2 = st.columns(2)
    lc1.link_button("View on NCBI Pathogen Detection", ncbi_url, use_container_width=True)
    lc2.link_button("Related activity on SalmonellaWatch-Ag", dashboard_url, use_container_width=True)


def render_serovar_section(s):
    if s.get("status") != "ok":
        st.warning(f"SISTR: {s.get('reason') or s.get('error', 'failed')}")
        return
    serovar = s.get("serovar") or "—"
    parts = serovar.split("|")
    primary = parts[0]
    rest = parts[1:] if len(parts) > 1 else []
    st.markdown(f"**Serovar (SISTR)**: {primary}")
    if rest:
        st.caption(f"Antigenically equivalent serovars (same formula): {', '.join(rest)}")
    af = s.get("antigenic_formula")
    if af:
        st.caption(f"Antigenic formula: {af}")
    if s.get("qc_status") and s["qc_status"] != "PASS":
        st.warning(f"SISTR QC: {s['qc_status']}")


def render_mlst_section(m):
    if m.get("status") != "ok":
        st.warning(f"MLST: {m.get('reason') or m.get('error', 'failed')}")
        return
    st.markdown(f"**Achtman MLST**: ST{m.get('st', '?')}")
    alleles = m.get("alleles") or []
    if alleles:
        st.caption("Alleles: " + ", ".join(alleles))


def render_plasmid_section(p):
    if p.get("status") != "ok":
        st.warning(f"PlasmidFinder: {p.get('reason') or p.get('error', 'failed')}")
        return
    hits = p.get("replicons") or []
    if not hits:
        st.markdown("**Plasmid replicons**: none detected")
        return
    st.markdown(f"**Plasmid replicons**: {len(hits)} detected")
    st.dataframe(hits, use_container_width=True, hide_index=True)


def render_pesi_section(pe):
    if pe.get("status") != "ok":
        st.warning(f"pESI: {pe.get('reason') or pe.get('error', 'failed')}")
        return
    verdict = pe.get("verdict", "—")
    n = pe.get("n_present", 0)
    panel = pe.get("n_markers_in_panel", "?")
    color_cls = "badge-ok" if "no pESI" in verdict.lower() or "unlikely" in verdict.lower() \
                else "badge-signal" if "possible" in verdict.lower() else "badge-soc"
    st.markdown(
        f"**pESI screen**: <span class='badge {color_cls}'>{verdict}</span> "
        f"({n} markers present of {panel} in panel)",
        unsafe_allow_html=True,
    )
    if pe.get("markers_present"):
        st.caption("Markers detected: " + ", ".join(pe["markers_present"]))
    st.markdown(
        '<div class="caveat-box">'
        "<strong>Manual review required.</strong> pESI is a ~280–300 kb plasmid that fragments "
        "across many short-read contigs. Any positive call here is a marker-supported inference, "
        "not a confirmed plasmid. Confirm by inspecting the plasmid contigs, ideally with a long-read "
        "or hybrid assembly."
        "</div>",
        unsafe_allow_html=True,
    )


def render_amrfinder_section(af):
    if af.get("status") != "ok":
        st.warning(f"AMRFinderPlus: {af.get('reason') or af.get('error', 'failed')}")
        return
    tabs = st.tabs([
        f"Virulence ({len(af.get('virulence', []))})",
        f"AMR ({len(af.get('amr', []))})",
        f"Stress ({len(af.get('stress', []))})",
        f"Point mutations ({len(af.get('point_mutations', []))})",
    ])
    for tab, key in zip(tabs, ["virulence", "amr", "stress", "point_mutations"]):
        with tab:
            rows = af.get(key) or []
            if not rows:
                st.caption("(none detected)")
                continue
            # show method column with color hint via a small legend
            st.dataframe(rows, use_container_width=True, hide_index=True)
    st.caption(
        "Method codes — **EXACTX** exact protein match · **BLASTX** translated nucleotide match · "
        "**PARTIALX** truncated (often at contig edge) · **INTERNAL_STOP** disrupted gene · "
        "**POINTX** resistance point mutation"
    )


def render_record(record, settings):
    sections = record.get("sections", {})
    st.markdown("---")
    st.subheader(record.get("assembly", "(unknown)"))

    # cluster card up top
    render_cluster_section(sections.get("cluster", {}), settings)

    # serovar / ST / plasmid in a row
    c1, c2 = st.columns(2)
    with c1:
        render_serovar_section(sections.get("sistr", {}))
        render_mlst_section(sections.get("mlst", {}))
    with c2:
        render_plasmid_section(sections.get("plasmid", {}))
        render_pesi_section(sections.get("pesi", {}))

    # AMRFinderPlus
    with st.expander("Virulence, AMR, stress, point mutations (AMRFinderPlus)", expanded=True):
        render_amrfinder_section(sections.get("amrfinder", {}))

    # raw JSON for the curious
    with st.expander("Raw JSON", expanded=False):
        st.json(record)


# ---------- UI ----------

st.title("🧬 SalmonellaWatch-Ag cluster matcher")
st.caption(
    "Identify the nearest *Salmonella enterica* SNP cluster and produce a full surveillance "
    "profile (serovar, ST, replicons, virulence, AMR, pESI screen) for one or more local "
    "assemblies. Files never leave your computer."
)

# sidebar: history + settings
with st.sidebar:
    st.subheader("Index status")
    build_date = S.get("index_build_date", "unknown")
    try:
        age = (datetime.now() - datetime.strptime(build_date, "%Y-%m-%d")).days
        color = "🟢" if age <= 90 else "🟡" if age <= 180 else "🔴"
        st.markdown(f"{color} Index built **{build_date}** ({age} days ago)")
    except Exception:
        st.markdown(f"Index built **{build_date}**")
    st.caption("Quarterly manual refresh recommended.")

    st.subheader("Recent runs")
    hist = list_history()
    if not hist:
        st.caption("No past runs yet.")
    else:
        for p in hist:
            label = f"{p.parent.name} · {p.stem}"
            if st.button(label, key=str(p), use_container_width=True):
                st.session_state["replay"] = str(p)

# main pane: file upload + run
uploaded = st.file_uploader(
    "Drop one or more .fna or .fasta files",
    type=["fna", "fasta", "fa"],
    accept_multiple_files=True,
)

run_clicked = st.button(
    f"Run on {len(uploaded)} file{'s' if len(uploaded) != 1 else ''}",
    type="primary",
    disabled=not uploaded,
)

# replay history if user clicked one
if "replay" in st.session_state:
    p = Path(st.session_state.pop("replay"))
    try:
        with open(p) as f:
            record = json.load(f)
        st.info(f"Showing saved result: {p.name}")
        render_record(record, S)
    except Exception as e:
        st.error(f"Could not load {p}: {e}")

# run pipeline
if run_clicked and uploaded:
    progress = st.progress(0.0)
    status = st.empty()
    results_container = st.container()

    tmp_dir = HERE / "_tmp_uploads"
    tmp_dir.mkdir(exist_ok=True)

    for i, uf in enumerate(uploaded, 1):
        status.markdown(f"**Processing {i}/{len(uploaded)}**: {uf.name}")
        tmp_path = tmp_dir / uf.name
        with open(tmp_path, "wb") as f:
            f.write(uf.getbuffer())

        record = run_panel(tmp_path, S["cpu"])
        if "error" in record:
            with results_container:
                st.markdown("---")
                st.subheader(uf.name)
                st.error(record["error"])
        else:
            save_result(uf.name, record)
            with results_container:
                render_record(record, S)

        try:
            tmp_path.unlink()
        except OSError:
            pass
        progress.progress(i / len(uploaded))

    status.success(f"Done — processed {len(uploaded)} assembl{'ies' if len(uploaded) != 1 else 'y'}.")
