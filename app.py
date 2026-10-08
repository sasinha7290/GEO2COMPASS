import concurrent.futures
import gc
import gzip
import io
import json
import os
import shutil
import re
import tarfile
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional
from urllib.parse import urljoin, parse_qs, urlparse
from urllib.request import Request, urlopen
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import uuid
from datetime import datetime, timezone

# PyArrow 25.0.0 can segfault when Streamlit initializes Arrow from a
# ScriptRunner thread. Use the system allocator even if the deployment
# environment does not define this variable .
os.environ.setdefault("ARROW_DEFAULT_MEMORY_POOL", "system")

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.types as pt
import pyarrow.dataset as ds
import duckdb
import GEOparse
import numpy as np
import pandas as pd
import requests
import pyarrow.parquet as pq
import streamlit as st
from bs4 import BeautifulSoup
from gprofiler import GProfiler

st.set_page_config(page_title="GEO-2-COMPASS")
st.title("GEO-2-COMPASS")

MAX_MATRIX_CELLS = 150_000_000
PREVIEW_ROWS = 1_000
PREVIEW_COLUMNS = 150

loaded_keys = [k for k in st.session_state.keys() if k.startswith("df_")]

with st.form("geo_accession_form"):
    requested_gse_id = st.text_input(
        "Enter GEO accession:",
        value=st.session_state.get("active_gse_id", ""),
        placeholder="e.g. GSE183620",
    )
    load_accession = st.form_submit_button("Load GEO metadata", type="primary")

if load_accession:
    requested_gse_id = requested_gse_id.strip().upper()
    if not re.fullmatch(r"GSE\d+", requested_gse_id):
        st.error("Enter a valid GEO Series accession, such as GSE183620.")
        st.stop()
    if requested_gse_id != st.session_state.get("active_gse_id"):
        for entry in st.session_state.get("result_lists") or []:
            for key in ("path", "modified_path"):
                if entry.get(key):
                    Path(entry[key]).unlink(missing_ok=True)
        st.session_state.clear()
        st.session_state["active_gse_id"] = requested_gse_id
        st.rerun()

gse_id = st.session_state.get("active_gse_id", "")
if not gse_id:
    st.info("Enter a GEO Series accession and select **Load GEO metadata**.")
    st.stop()

if "run_pipeline" not in st.session_state:
    st.session_state.run_pipeline = False
if "selected_columns" not in st.session_state:
    st.session_state.selected_columns = []
if "selected_gpl" not in st.session_state:
    st.session_state.selected_gpl = []
if "gpl_data" not in st.session_state:
    st.session_state.gpl_data = {}
if "current_gpl" not in st.session_state:
    st.session_state.current_gpl = {}
if "dp_text" not in st.session_state:
    st.session_state.dp_text = ""
if "char_df" not in st.session_state:
    st.session_state.char_df = pd.DataFrame()
if "record_log" not in st.session_state:
    st.session_state.record_log = []

def log_record(category, message, **details):
    if "record_log" not in st.session_state:
        st.session_state.record_log = []
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "category": category,
        "message": message,
    }
    if details:
        entry["details"] = details
    st.session_state.record_log.append(entry)

def fetch_metadata(gse_id: str, gse) -> dict:
    try:
        gpl_ids = list(gse.gpls.keys())
        gpl_titles = [
            gse.gpls[gpl_id].metadata.get('title', [''])[0]
            for gpl_id in gpl_ids
        ]

        gsm_list = list(gse.gsms.keys())

        gsm_to_gpl = {}
        gsm_gpl_dict = {gpl_id: [] for gpl_id in gpl_ids}

        gsm_gpl_dict["Unknown"] = []

        for gsm_id, gsm_obj in gse.gsms.items():
            platform_list = gsm_obj.metadata.get('platform_id', [])
            associated_gpl = platform_list[0] if platform_list else ""
            gsm_to_gpl[gsm_id] = associated_gpl

            if associated_gpl in gsm_gpl_dict:
                gsm_gpl_dict[associated_gpl].append(gsm_id)
            elif associated_gpl:
                gsm_gpl_dict[associated_gpl] = [gsm_id]
            else:
                gsm_gpl_dict["Unknown"].append(gsm_id)

        gse_types = gse.metadata.get('type', [])
        gse_types_str = " ".join(gse_types).lower()
        study_type = "Microarray" if "array" in gse_types_str else "RNA-seq"

        taxon = ""
        if 'organism' in gse.metadata and gse.metadata['organism']:
            taxon = gse.metadata['organism'][0]
        elif gse.gsms:
            first_gsm = next(iter(gse.gsms.values()))
            organism_list = first_gsm.metadata.get('organism_ch1', [])
            if organism_list:
                taxon = organism_list[0]

        return {
            "title":      gse.metadata.get('title', [gse_id])[0],
            "gse_id": gse_id,
            "type":       study_type,
            "gsm_ids":    gsm_list,
            "gpl_ids":    gpl_ids,
            "gpl_titles": gpl_titles,
            "gsm_gpl_dict": gsm_gpl_dict,
            "taxon":      taxon,
            "error":      None,
        }

    except Exception as e:
        return {"error": str(e)}

@st.cache_resource(
    show_spinner="Downloading and parsing GEO metadata…",
    ttl=300,
    max_entries=2,
)
def get_gse(gse_id):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/soft/{gse_id}_family.soft.gz"
    headers = {
        "User-Agent": "Mozilla/5.0 (compatible; GEO-2-COMPASS/1.0; contact: your_email@example.com)"
    }
    print(f" [Downloading] Fetching data directly from: {url}")

    try:
        with requests.get(
            url,
            stream=True,
            headers=headers,
            timeout=(30, 300),
        ) as response:
            response.raise_for_status()

            with tempfile.NamedTemporaryFile(
                suffix=".soft.gz",
                delete=False,  ) as temp_file:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        temp_file.write(chunk)
                temp_file.flush()

            try:
                print("Parsing data with GEOparse from temporary file...")
                return GEOparse.get_GEO(
                    filepath=temp_file.name,
                    geotype="GSE",
                )
            finally:
                os.unlink(temp_file.name)
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Download failed: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Failed to parse GEO data: {e}") from e

def get_dp_and_char(gse, meta):
    first_gsm = gse.gsms[meta["gsm_ids"][0]]
    dp_text = first_gsm.metadata.get("data_processing", [""])[0]
    st.session_state["dp_text"] = dp_text
    char_list = set()
    for gsm in meta["gsm_ids"]:
        gsm_data = gse.gsms[gsm]
        for x in gsm_data.metadata.get("characteristics_ch1", []):
            if ": " in x:
                char_list.add(x.split(": ", 1)[0])
            else:
                char_list.add(x)

    char_list = list(char_list)

    char_df = pd.DataFrame(columns=char_list)
    for gsm in meta["gsm_ids"]:
        gsm_data = gse.gsms[gsm]
        new_row = {}

        for x in gsm_data.metadata.get("characteristics_ch1", []):
            if ": " in x:
                key, value = x.split(": ", 1)
                new_row[key] = value
            else:
                new_row[x] = x
        if "characteristics_ch1" in gsm_data.metadata:
            char_df.loc[gsm_data.metadata["geo_accession"][0]] = new_row
        for key in ["treatment_protocol_ch1", "growth_protocol_ch1", "organism_ch1", "source_name_ch1", "title"]:
            if key in gsm_data.metadata:
                char_df.loc[gsm_data.metadata["geo_accession"][0], key] = gsm_data.metadata[key][0]
    st.session_state["char_df"] = char_df.astype("string")

try:
    GLOBAL_GSE = get_gse(gse_id)
except RuntimeError as exc:
    st.error(str(exc))
    st.stop()

meta = fetch_metadata(gse_id, GLOBAL_GSE)

if meta.get("error"):
    st.error(f"Failed to load metadata: {meta['error']}")
    st.stop()

if not meta.get("gsm_ids") or not meta.get("gpl_ids"):
    st.error("The GEO record does not contain usable sample or platform metadata.")
    st.stop()


st.subheader(meta["title"])
st.caption(
    f"Type: **{meta['type']}** · "
    f"Samples: **{len(meta['gsm_ids'])}** · "
    f"Taxon: {meta['taxon']}"
)

st.session_state.selected_gpl = meta["gpl_ids"]

get_dp_and_char(GLOBAL_GSE, meta)

# # # # # # # # # # # # # # # #
# NORMALIZATION DETECTION # # #
# # # # # # # # # # # # # # # #

@dataclass
class NormInfo:
    norm_type: str
    is_log: bool
    matched_term: str
    source: str

def _normalize_for_match(text):
    return re.sub(r"[^a-zA-Z0-9]+", " ", text or "")

_NORM_RULES: list[tuple[str, "re.Pattern", bool]] = [
    ("median_of_ratios", re.compile(r"\bdeseq2?\b|\bmedian of ratios\b|\bsize factors?\b", re.I), False),
    ("tmm",              re.compile(r"\btmm\b|\bedger\b|\btrimmed mean\b", re.I), False),
    ("upper_quartile",   re.compile(r"\bupper quartile\b|\buq\b", re.I), False),
    ("voom",             re.compile(r"\bvoom\b", re.I), True),
    ("vst",              re.compile(r"\bvst\b|\bvariance stabili[sz]", re.I), True),
    ("rlog",             re.compile(r"\brlog\b|\bregulari[sz]ed log\b", re.I), True),
    ("tpm",              re.compile(r"\btpm\b|\btranscripts? per million\b", re.I), False),
    ("rpkm_fpkm",        re.compile(r"\brpkm\b|\bfpkm\b|\b(reads|fragments) per kilobase\b", re.I), False),
    ("cpm",              re.compile(r"\bcpm\b|\bcounts? per million\b", re.I), False),
    ("quantile",         re.compile(r"\bquantile\b", re.I), False),
    ("rma",              re.compile(r"\brma\b|\bfrma\b|\bgcrma\b|\brobust multi ?array average\b", re.I), True),
    ("mas5",             re.compile(r"\bmas ?5\b", re.I), False),
    ("plier",            re.compile(r"\bplier\b", re.I), True),
    ("loess",            re.compile(r"\bloess\b|\blowess\b", re.I), False),
    ("beadstudio",       re.compile(r"\bbeadstudio\b|\bgenomestudio\b", re.I), True),
    ("geometric",        re.compile(r"\bgeometric\b", re.I), False),
    ("raw_counts",       re.compile(r"\braw\b|\bunnormali[sz]ed\b|\bunnorm\b|\bhtseq ?count\b|\bfeaturecounts\b|\brsem\b|\bexpected count\b", re.I), False),
]

_EXPLICIT_LOG_RE = re.compile(r"\blog ?2\b|\blog ?10\b|\blog transform|\blog scal", re.I)
_EXPLICIT_NOLOG_RE = re.compile(r"\bnot log|\bnon ?log\b|\blinear scale\b|\braw scale\b|\buntransformed\b", re.I)

def _scan_norm_keywords(raw_text: str):
    text = _normalize_for_match(raw_text)
    best = None
    for norm_type, pattern, default_log in _NORM_RULES:
        for m in pattern.finditer(text):
            if best is None or m.start() > best[3]:
                best = (norm_type, m.group(0), default_log, m.start())
    return best

def detect_normalization(dp_text: str = "", filename: str = "") -> NormInfo:
    hit = _scan_norm_keywords(dp_text)
    source = "data_processing_text"
    if hit is None and filename:
        hit = _scan_norm_keywords(filename)
        source = "filename"
    if hit is None:
        return NormInfo("unknown", False, "", "default")

    norm_type, matched_term, default_log, _ = hit
    combined = _normalize_for_match(f"{dp_text} {filename}")
    if _EXPLICIT_NOLOG_RE.search(combined):
        is_log = False
    elif _EXPLICIT_LOG_RE.search(combined):
        is_log = True
    else:
        is_log = default_log

    return NormInfo(norm_type, is_log, matched_term, source)


def to_cpm(df):
    numeric_df = df.apply(pd.to_numeric, downcast="float", errors="coerce")
    col_sums = numeric_df.sum(axis=0).replace(0, np.nan)
    return numeric_df.divide(col_sums, axis=1) * 1e6

def reduce_matrix_memory(df: pd.DataFrame) -> pd.DataFrame:
    """Downcast numeric columns so one large matrix fits a small web worker."""
    for column in df.select_dtypes(include=[np.number]).columns:
        if pd.api.types.is_integer_dtype(df[column]):
            df[column] = pd.to_numeric(df[column], downcast="integer")
        else:
            df[column] = pd.to_numeric(df[column], downcast="float")
    return df

# # # # # # # # # # # # # # # #
# GENE SYMBOL ANNOTATION  # # #
# # # # # # # # # # # # # # # #

_REJECT_COMBINED = re.compile(
    r"^(?:NM|NR|NP|XM|XR|XP|NG|NT|NC)_\d+(?:\.\d+)?$|"
    r"^AFFX-.*_at$|"
    r"^.*_at$|"
    r"^[A-Z]{1,2}\d{6,}(?:\.\d+)?$|"
    r"^GenMAPP|"
    r"^ILMN_\d+$|"
    r"^(?:ENSG|ENST|ENSMUSG|ENSMUST|ENSP)\d+|"
    r"^\d+_(?:[a-z]_at|at|x_at|s_at)$|"
    r"^A_\d+_P\d+$|"
    r"^chr|"
    r"^\d+[pq]\d+|"
    r"^GO:\d+$|"
    r"\.\d+$|"
    r"^\d|\s|\.|_\d{4,}|scl\d+\.|RefSeq|\+|"
    r"^[agctAGCT]*$|"
    r"^.{0,3}$",
    re.IGNORECASE,
)

_ACCEPT_COMBINED = re.compile(r"^\d{7}[A-Za-z]\d{2}[Rr][Ii][Kk]$")


def is_gene_symbol(value):
    s_val = str(value)
    if "//" in s_val:
        parts = s_val.split("//")
        v = parts[1].strip() if len(parts) > 1 else parts[0].strip()
    else:
        v = s_val.strip()

    if _REJECT_COMBINED.search(v):
        return bool(_ACCEPT_COMBINED.search(v))

    return True

def too_homogenous(col):
    if len(col) <= 0:
        return True
    if "//" in str(list(col)[0]):
        col = list(set([
            parts[1]
            for c in col
            if len(parts := re.split(r'///|//', str(c), maxsplit=1)) > 1
        ]))
        if len(col) < 100:
            return True

    results = []
    for c in col:
        result = re.split(r'(?<=\D)(?=\d)', str(c))
        results.append(result[0])

    if len(list(set(results))) < 30:
        return True

    col = list(set(col))
    return len(col) < 100

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

def fetch_hgnc_table() -> pd.DataFrame:
    df = pd.read_csv(HGNC_URL, sep="\t", dtype=str, low_memory=False)
    return df[["hgnc_id", "symbol", "prev_symbol", "alias_symbol"]]

def gene_convert(gpl_df, gse, best_col, symbol_col):
    for gpl_name, gpl in getattr(gse, "gpls", {}).items():
        try:
            species = gpl.metadata.get('organism', [])[0]
        except (AttributeError, TypeError, KeyError, IndexError):
            return gpl_df, symbol_col


    if (any(s in species.lower() for s in ["homo","sapiens"])):
        species = "hsapiens"
    elif (any(s in species.lower() for s in ["musculus","mus"])):
        species = "mmusculus"
    else:
        return gpl_df, symbol_col

    target_namespace = "HGNC" if species == "hsapiens" else "MGI"

    id = str(gpl_df[best_col].iloc[0])
    numeric_namespace = None
    if (bool(re.match(r"^\d+(_at)?$", id))):
        numeric_namespace = "ENTREZGENE_ACC"
    if (bool(re.match(r"^ENSG\d{11}(\.\d+)?$", id)) or bool(re.match(r"^ENSMUSG\d{11}(\.\d+)?$", id))):
        numeric_namespace = "ENSEMBL"
    if (id.strip()[:3].lower() == "eg:"):
        numeric_namespace = "ENTREZGENE_ACC"


    if numeric_namespace == None:
        return gpl_df, symbol_col
    else:
        gpl_df["new_id_col"] = (
            gpl_df[best_col].astype(str)
            .str.replace(r"_at$", "", regex=True)
            .str.replace(r"\.\d+$", "", regex=True)
        )
        if (str(gpl_df[best_col].iloc[0]).strip()[:3].lower() == "eg:"):
            gpl_df["new_id_col"] = gpl_df[best_col].astype(str).str.split(':').str[1]

        try:
            gp = GProfiler(return_dataframe=True)
            results = gp.convert(organism=species, query=list(gpl_df["new_id_col"]),
                                numeric_namespace=numeric_namespace, target_namespace=target_namespace)
            if results.empty or "name" not in results.columns:
                log_record(
                    "gene_annotation", "g:Profiler conversion returned no results; original IDs kept.",
                    species=species, numeric_namespace=numeric_namespace, target_namespace=target_namespace,
                )
                return gpl_df, symbol_col
        except Exception as e:
            log_record(
                "gene_annotation", "g:Profiler conversion failed; original IDs kept.",
                error=str(e), species=species,
            )
            return gpl_df, symbol_col

        results_deduped = results.drop_duplicates(subset="incoming", keep="first")
        mapping = results_deduped.set_index("incoming")["name"]

        gpl_df["gene_symbol"] = [mapping.get(id_, float("nan")) for id_ in gpl_df["new_id_col"]]

        mapped_count = int(gpl_df["gene_symbol"].notna().sum())
        log_record(
            "gene_annotation", "Gene symbols assigned via g:Profiler conversion.",
            species=species, numeric_namespace=numeric_namespace, target_namespace=target_namespace,
            mapped=mapped_count, total=len(gpl_df),
            mapping_rate=round(mapped_count / len(gpl_df), 3) if len(gpl_df) else 0,
        )

        return gpl_df, "gene_symbol"

def select_gpl(results):
    gene_columns = []
    for x in results.values():
        mapped = x[3] if len(x) > 3 else []
        if mapped:
            gene_columns.append(mapped)
        else:
            gpl_df, symbol_col = x[0], x[2]
            if symbol_col is not None and symbol_col in getattr(gpl_df, "columns", []):
                gene_columns.append(gpl_df[symbol_col].dropna().astype(str).tolist())
            else:
                gene_columns.append([])

    intersection = set()
    union = set()


    for g in gene_columns:
        if len(intersection) == 0:
            intersection = set(g)
            union = union | set(g)
        else:
            intersection = intersection & set(g)
            union = union | set(g)
    if len(union) == 0:
        return results

    print(len(intersection) / len(union))

    if len(intersection) / len(union) > 0.65:
        log_record(
            "gpl_selection",
            "Platforms merged automatically (gene sets sufficiently overlapping).",
            overlap_ratio=round(len(intersection) / len(union), 3),
            gpl_ids=list(results.keys()),
        )
        return results
    else:
        if len(meta["gpl_ids"]) > 1:
            selected_gpl = st.selectbox(
                "Multiple platforms detected. Which platform would you like to use?",
                tuple(f"{meta['gpl_ids'][i]}:  {meta['gpl_titles'][i]}" for i in range(len(meta['gpl_ids']))),
                index=0,
                placeholder="Select GPL"
            )
            if selected_gpl:
                for key in list(st.session_state.keys()):
                    if str(key).startswith(("survival")):
                        del st.session_state[key]
                gc.collect()
                selected_gpl = selected_gpl.split(":  ")[0]
                meta['gsm_ids'] = meta["gsm_gpl_dict"][selected_gpl]
                meta['gpl_ids'] = [selected_gpl]
                log_record(
                    "gpl_selection",
                    "User selected a single platform; platforms were not merged.",
                    overlap_ratio=round(len(intersection) / len(union), 3),
                    selected_gpl=selected_gpl,
                    available_gpls=meta["gpl_ids"],
                )
            else:
                st.info("Select one platform before building the matrix.")
                st.stop()
        else:
            selected_gpl = meta["gpl_ids"][0]
        return {selected_gpl: results[selected_gpl]}


def get_gene_symbol_column(counts_df, all_selected_gpls, probe_ids=None):
    gse = GLOBAL_GSE
    results = {}

    for selected_gpl in all_selected_gpls:
        gpl_df = None

        for gpl_name, gpl in getattr(gse, "gpls", {}).items():
            if gpl_name == selected_gpl:
                gpl_df = gpl.table

        if gpl_df is None:
            raise ValueError(f"Platform '{selected_gpl}' not found in this GEO series; skipping gene symbol annotation.")

        if probe_ids is not None:
            index_set = {_normalize_id(v) for v in probe_ids}
        else:
            index_set = {_normalize_id(v) for v in counts_df.index}

        if gpl_df is None or gpl_df.empty:
            gpl_df = pd.DataFrame({"id": list(index_set)})
            gpl_df["gene_symbol"] = gpl_df["id"]

            gpl_df, symbol_col = gene_convert(gpl_df, gse, "id", "gene_symbol")
            results[selected_gpl] = [gpl_df, "id", symbol_col]
            continue

        best_col, best_score = _find_best_id_col(index_set, gpl_df)

        if best_score < 0.01:
            print(f"  [warn] best overlap is only {best_score:.3f} — index may not match any GPL column")

        if gpl_df is not None:
            symbol_col = None
            symbol_col_score = 0
            for col in gpl_df.columns:
                values = gpl_df[col].replace("---", pd.NA).dropna()
                if too_homogenous(values):
                    continue
                col_score = values.apply(is_gene_symbol).mean()
                print(f"  symbol candidate '{col}': {col_score:.3f}")
                if col_score > symbol_col_score:
                    symbol_col_score = col_score
                    symbol_col = col
            if symbol_col_score < 0.30: #arbitrary
                gpl_df, symbol_col = gene_convert(gpl_df, gse, best_col, symbol_col)

            results[selected_gpl] = [gpl_df, best_col, symbol_col]
        else:
            raise ValueError(f"No platform data found for {gse}")

    return results


def _strip_version(s: str) -> str:
    return re.sub(r'\.\d+$', '', s)

def _normalize_id(s: str) -> str:
    s = str(s).strip().upper()
    if re.match(r'^\d+\.0$', s):
        s = s[:-2]
    return s

def _overlap_score(index_set: set[str], col_vals: pd.Series) -> float:
    col_set_raw = set(col_vals.apply(_normalize_id))

    # exact (normalized)
    exact = len(index_set & col_set_raw) / len(index_set)
    if exact > 0:
        return exact

    # version-stripped
    index_stripped = {_strip_version(v) for v in index_set}
    col_stripped   = {_strip_version(v) for v in col_set_raw}
    return len(index_stripped & col_stripped) / len(index_set)


def _find_best_id_col(index_set: set[str], gpl_df: pd.DataFrame):
    best_col   = None
    best_score = -1

    for col in gpl_df.columns:
        score = _overlap_score(index_set, gpl_df[col])
        print(f"  id candidate '{col}': {score:.3f}")
        if score > best_score:
            best_score = score
            best_col   = col

    print(f"  → best id col: '{best_col}' (score {best_score:.3f})")
    return best_col, best_score

# # # # # # # # # # # # # # # #
# RNASEQ DATASET HANDLING # # #
# # # # # # # # # # # # # # # #

@dataclass
class FileMeta:
    """Metadata for one downloadable supplementary file."""
    url:      str
    filename: str
    is_tar:   bool = False
    ncbi_data: bool = False
    normalization: str = "unknown"
    is_log: bool = False

    @property
    def is_tabular(self) -> bool:
        return _is_tabular(self.filename)

@dataclass
class GseResult:
    """Everything produced for one GEO accession."""
    accession: str
    normalization_type: str
    effective_norm: str

    counts_df:   Optional[pd.DataFrame] = field(default=None, repr=False)
    norm_df:     Optional[pd.DataFrame] = field(default = None, repr = False)

    def save(self, output_dir: str | Path = ".") -> list[Path]:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        written: list[Path] = []

        target_df   = self.norm_df if self.norm_df is not None else self.counts_df
        target_name = "log2CPM" if self.norm_df is not None else "counts"

        if target_df is not None:
            p = output_dir / f"{self.accession}_{target_name}.txt"
            target_df.to_csv(p, sep="\t")
            print(f"✓ {target_name} matrix → {p}  "
                  f"({target_df.shape[0]} genes × {target_df.shape[1]} samples)")
            written.append(p)

        mp = output_dir / f"{self.accession}_meta.json"
        with open(mp, "w") as f:
            json.dump({
                "accession":          self.accession,
                "normalization_type": self.normalization_type,
                "effective_norm":     self.effective_norm,
            }, f, indent=2)
        written.append(mp)
        return written

    def __repr__(self) -> str:
        df = self.norm_df if self.norm_df is not None else self.counts_df
        shape = f"{df.shape[0]}g x {df.shape[1]}s" if df is not None else "no matrix"
        return (f"GseResult(acc={self.accession!r}, "
                f"orig_norm={self.normalization_type!r}, "
                f"effective={self.effective_norm!r}, ")

def _is_gzip(path: Path) -> bool:
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"

def _gunzip_to_disk(path: Path) -> Path:
    """Stream-decompress without holding the file in memory."""
    out = path.with_suffix(".plain")
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst, 1024 * 1024)
    return out

def _peek(path: Path, n: int = 16384) -> bytes:
    with open(path, "rb") as f:
        return f.read(n)

MAX_DOWNLOAD_BYTES = 5 * 1024**3
DOWNLOAD_DIR = Path(tempfile.gettempdir()) / "geo_downloads"

def _download_bytes(url: str, cap=MAX_DOWNLOAD_BYTES) -> Optional[Path]:
    """Stream url to a temp file. Returns its Path, or None if over the cap."""
    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    dest = DOWNLOAD_DIR / f"{uuid.uuid4().hex}.download"
    written = 0
    try:
        with requests.get(url, stream=True, timeout=(30, 300)) as r:
            r.raise_for_status()
            declared = int(r.headers.get("Content-Length") or 0)
            if declared > cap:
                return None
            with open(dest, "wb") as f:
                for chunk in r.iter_content(1024 * 1024):
                    written += len(chunk)
                    if written > cap:
                        raise _TooBig()
                    f.write(chunk)
        return dest
    except _TooBig:
        dest.unlink(missing_ok=True)
        return None
    except Exception:
        dest.unlink(missing_ok=True)
        raise

class _TooBig(Exception):
    pass


def _download_all(urls: list[str]) -> dict[str, Optional[Path]]:
    results: dict[str, Optional[Path]] = {}
    if not urls:
        return results

    workers = min(1, len(urls))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        fut_to_url = {ex.submit(_download_bytes, u): u for u in urls}
        for fut in concurrent.futures.as_completed(fut_to_url):
            url = fut_to_url[fut]
            fname = url.split("/")[-1]
            try:
                res = fut.result()
                results[url] = res
                if res is not None:
                    print(f"    ✓ {fname} ({res.stat().st_size / 1e6:.1f} MB on disk)")
                else:
                    print(f"    ✗ {fname}: Failed or exceeded size limit")
            except Exception as e:
                print(f"    ✗ {fname}: {e}")
    return results

#simple helper functions
TABULAR_EXTS = {".txt", ".csv", ".tsv", ".tab", ".xlsx", ".xls"}
def _is_tabular(name: str) -> bool:
    lower = name.lower()
    return any(lower.endswith(e) or lower.endswith(e + ".gz") for e in TABULAR_EXTS)


def _decompress(raw: bytes, name: str) -> tuple[bytes, str]:
    if name.lower().endswith(".gz"):
        return gzip.decompress(raw), name[:-3]
    return raw, name


def detect_delimiter(
    raw_bytes: bytes,
    candidates=(",", "\t", ";", "|", "\n", "\r\n", "\r", " ", "~"),
    sample_size: int = 16384,
):
    """Detects both column separator and record/row terminator from raw bytes.

    Evaluates candidate pairs (col_sep, row_sep) by scoring grid uniformity.
    """
    text = raw_bytes[:sample_size].decode("utf-8", errors="ignore")

    best_pair = (",", "\n")
    best_score = -1.0

    # Evaluate every combination of column separator S and record separator R
    for r_sep in candidates:
        # Split records using candidate row separator
        records = [r for r in text.split(r_sep) if r.strip()][:30]
        if len(records) < 2:
            continue

        for c_sep in candidates:
            if c_sep == r_sep:
                continue

            # Count fields per record
            row_lengths = [len(rec.split(c_sep)) for rec in records]
            if not row_lengths:
                continue

            avg_cols = sum(row_lengths) / len(row_lengths)

            # Ignore single-column results or empty splits
            if avg_cols <= 1.05:
                continue

            # Calculate variance in column count across records
            variance = sum((x - avg_cols) ** 2 for x in row_lengths) / len(
                row_lengths
            )

            # Score prioritizes high field count with perfect uniformity (0 variance)
            score = (avg_cols * len(records)) / (1.0 + (variance * 10.0))

            if score > best_score:
                best_score = score
                best_pair = (c_sep, r_sep)

    return best_pair


#tested
def _scrape_geo_download_page(accession):
    page_url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession}"
    print(f"  Scraping GEO download page: {page_url}")

    try:
        r = requests.get(page_url, timeout=60)
        r.raise_for_status()
    except requests.exceptions.Timeout:
        st.error("GEO download page timed out. Try again in a moment.")
        st.stop()
    except requests.exceptions.RequestException as e:
        st.error(f"Could not reach GEO download page: {e}")
        st.stop()

    soup = BeautifulSoup(r.text, "html.parser")

    seen = set()
    files = []

    for a in soup.find_all("a", href=True):
        href = a["href"]

        if "file=" in href:
            filename = href.split("file=")[-1].split("&")[0]
        elif "suppl/" in href:
            filename = href.split("suppl/")[-1]
        else:
            continue

        if not filename or filename in seen:
            continue
        seen.add(filename)

        full_url = urljoin(page_url, href)
        lower    = filename.lower()
        is_tar   = lower.endswith(".tar") or lower.endswith(".tar.gz")
        ncbi_data = "file=" in href

        norm_info = detect_normalization(filename=filename)
        normalization = norm_info.norm_type

        print(f"\n\n\n {filename} \t {normalization} \n\n\n")


        if not (is_tar or _is_tabular(filename)) or "annot" in filename:
            print(f"    [skip] {filename}")
            continue

        files.append(FileMeta(
            url=full_url, filename=filename, is_tar=is_tar, ncbi_data=ncbi_data,
            normalization=normalization, is_log=norm_info.is_log,
        ))
        print(f"    found: {filename}")

    print(f"{len(files)} candidates on download page")
    return files

#constructing the dataframe

def _parse_tabular_series(raw: bytes, filename: str) -> tuple[pd.Series, str]:
    sep = detect_delimiter(raw)[0]
    df  = pd.read_csv(
        io.StringIO(raw.decode("utf-8", errors="replace")),
        sep=sep, comment="#", header=0,
    )

    mask = df.iloc[:, 0].astype(str).str.startswith("__")
    df = df[~mask.astype(bool)]

    num_cols = [c for c in df.columns if pd.to_numeric(df[c], errors="coerce").notna().all()]
    str_cols = [c for c in df.columns if c not in num_cols]

    # Gene index from first non-numeric column
    gene_index = df[str_cols[0]].astype(str) if str_cols else df.index.astype(str)

    if len(num_cols) == 0:
        raise ValueError(f"No numeric columns in {filename}")
    else:
        count_col = num_cols[0]

    s = pd.to_numeric(
        df[count_col],
        errors="coerce",
        downcast="integer",
    )
    s.index = gene_index

    stem = re.sub(r"\.gz$", "", filename, flags=re.IGNORECASE)
    stem = re.sub(r"\.(csv|tsv|txt|tab)$", "", stem, flags=re.IGNORECASE)
    return s, stem


def _build_matrix_from_tar(url_bytes, file_meta):
    series_list = []

    for meta in file_meta:
        path = url_bytes.get(meta.url)
        if path is None:
            continue
        with tarfile.open(name=path, mode="r:*") as tf:
            for member in tf:
                if not member.isfile():
                    continue
                mname = Path(member.name).name

                if not _is_tabular(mname):
                    continue
                extracted = tf.extractfile(member)
                if extracted is None:
                    continue
                try:
                    mbytes, bare = _decompress(extracted.read(), mname)
                    s, sname = _parse_tabular_series(mbytes, bare)
                    s.name = sname
                    series_list.append(s)
                    del mbytes
                except Exception as e:
                    print(f"    [warn] {mname}: {e}")
                    return pd.DataFrame()

    if not series_list:
        return pd.DataFrame()

    print(f"\n  Merging {len(series_list)} sample series from TAR…")
    df = pd.concat(series_list, axis=1, join="outer")
    df = df.apply(lambda c: pd.to_numeric(c, errors="coerce", downcast="float"))
    df.index.name = "gene_id"
    print(f"  ✓ Matrix: {df.shape[0]} genes × {df.shape[1]} samples")
    return df

def _build_matrix_from_flat_files(url_bytes, file_meta):
    dfs = []

    for meta in file_meta:
        path = url_bytes.get(meta.url)
        if path is None:
            continue
        filename = meta.filename
        plain = path
        try:
            if _is_gzip(path):
                plain = _gunzip_to_disk(path)
                if filename.endswith(".gz"):
                    filename = filename[:-3]

            lower = filename.lower()
            if lower.endswith(".xls") or lower.endswith(".xlsx"):
                df = pd.read_excel(plain, index_col=0)
            else:
                sr = detect_delimiter(_peek(plain))
                df = pd.read_csv(
                    plain, sep=sr[0],
                    lineterminator=sr[1],
                    index_col=0, on_bad_lines="skip",
                )
                numeric_cols = df.select_dtypes(include=[np.number]).columns
                df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, downcast="float")


            if df.columns.str.startswith('Unnamed').all():
                df.columns = df.iloc[0]

            df = df.apply(pd.to_numeric, downcast = "float", errors = "coerce")
            df.index.name = "gene_id"

            df = df[~df.index.astype(str).str.startswith("__")]

            frac_numeric = df.notna().mean()
            df = df.loc[:, frac_numeric > 0.3]



            dfs.append(df)

        except Exception as e:
            print(f"Error: {e}")
            return pd.DataFrame()
        finally:
            if plain != path:
                plain.unlink(missing_ok=True)

    if not dfs:
        return pd.DataFrame()

    if len(dfs) == 1:
        return dfs[0]

    print("Merging...")
    merged = pd.concat(dfs, axis=1, join = "outer")
    return merged


def _fetch_counts_df(accession, file_meta, url_bytes):

    tar_files  = [m for m in file_meta if m.is_tar]
    flat_files = [m for m in file_meta if not m.is_tar]

    if tar_files:
        return _build_matrix_from_tar(url_bytes, tar_files)
    elif flat_files:
        return _build_matrix_from_flat_files(url_bytes, flat_files)
    else:
        return pd.DataFrame()

def _annotate_counts(accession, selected, counts_df, all_selected_gpls):
    if selected[0].ncbi_data:
        raw = b""
        with gzip.open("annot/Human.GRCh38.p13.annot.tsv.gz", "rb") as f:
            raw = f.read()

        annot_df = pd.read_csv(
            io.BytesIO(raw), sep= "\t", on_bad_lines="skip",
        )

        id_col = "GeneID"
        symbol_col = "Symbol"
        results = {"ncbi": [annot_df, id_col, symbol_col]}
    else:
        try:
            results = get_gene_symbol_column(counts_df, all_selected_gpls)
        except ValueError as e:
            return counts_df, {}

    mapping = {}
    orig_index_str = counts_df.index.astype(str)
    orig_index_lower = orig_index_str.str.lower()

    for gpl, item in results.items():
        annot_df, id_col, symbol_col = item[0], item[1], item[2]

        gpl_map = dict(zip(
            annot_df[id_col].astype(str).str.lower(),
            annot_df[symbol_col].astype(str)
        ))

        mapping.update(gpl_map)

        mapped_genes = orig_index_lower.map(gpl_map).dropna().tolist()

        if isinstance(item, tuple):
            results[gpl] = list(item) + [mapped_genes]
        else:
            results[gpl].append(mapped_genes)

    new_index = orig_index_lower.map(mapping)
    counts_df.index = new_index.where(
        new_index.notna() & (new_index != "nan"), orig_index_str
    )
    counts_df.index.name = "gene_id"
    print(f"  ✓ Index remapped across {len(results)} GPL platform(s)")

    return counts_df, results


def _read_dp_text(geo, accession):
    if accession.startswith("GSE"):
        gsm_dict = geo.gsms
    else:
        gsm_dict = {accession: geo}

    first_gsm = next(iter(gsm_dict.values()))
    return "\n".join(first_gsm.metadata.get("data_processing", []))

def fetch_and_normalize(
    accession:     str,
    gsm_ids:      Optional[list[str]] = None,
    geo_cache_dir: str | Path = "./geo_cache",
    save_output:   bool       = False,
    output_dir:    str | Path = "./geo_output",
) -> list[dict[str, Any]]:

    geo_cache_dir = Path(geo_cache_dir)
    geo_cache_dir.mkdir(parents=True, exist_ok=True)
    accession = accession.strip().upper()

    print(f"\n{'═'*60}")
    print(f"  GEO accession: {accession}")
    print(f"{'═'*60}")
    dp_text  = st.session_state["dp_text"]

    candidates = _scrape_geo_download_page(accession)
    dp_norm_info = detect_normalization(dp_text=dp_text)
    if not candidates:
        return []

    results_meta = []
    merged_gpl_data = {}

    all_selected_gpls = st.session_state.selected_gpl

    def process_candidate(c):
        selected = [c]
        if c.normalization == "unknown":
            c.normalization = dp_norm_info.norm_type
            c.is_log = dp_norm_info.is_log

        log_record(
            "normalization", "Normalization detected from source metadata.",
            norm_type=c.normalization, is_log=c.is_log
        )

        print(f"\n  Processing: {c.filename}")

        url_bytes = _download_all([c.url])
        dl_path = url_bytes.get(c.url)

        if dl_path is None:
            print(f"    Exceeded file size, skipping {c.url}")
            return None

        try:
            counts_df = _fetch_counts_df(accession, selected, url_bytes)
        finally:
            dl_path.unlink(missing_ok=True)

        gc.collect()

        if counts_df.empty:
            return None

        if counts_df.size > MAX_MATRIX_CELLS:
            print(f"    Skipping {c.filename}: {counts_df.size:,} cells exceeds "
                  f"{MAX_MATRIX_CELLS:,}-cell limit")
            return None

        counts_df, gpl_results = _annotate_counts(accession, selected, counts_df, all_selected_gpls)

        retained_cols = [
            col for col in counts_df.columns
            if "gsm" not in col.lower() or any(gsm in col for gsm in gsm_ids)
        ]

        counts_df = counts_df[retained_cols]
        counts_df.index.name = None
        counts_df.insert(0, "Name", counts_df.index)
        counts_df.reset_index(drop=True, inplace=True)

        cache_id = uuid.uuid4().hex
        cache_path = geo_cache_dir / f"{accession}_{cache_id}.parquet"
        modified_path = geo_cache_dir / f"{accession}_{cache_id}_modified.parquet"

        # buffer = io.BytesIO()
        # counts_df.to_parquet(buffer, index=False)
        # buffer.seek(0)

        # cache.add(cache_path, buffer.getvalue())
        # cache.add(modified_path, buffer.getvalue())

        counts_df.to_parquet(cache_path)
        shutil.copy(cache_path, modified_path)

        meta_entry = {
            "path": cache_path,
            "modified_path": modified_path,
            "filename": c.filename,
            "normalization_type": c.normalization,
            "is_log": c.is_log,
            "n_genes": counts_df.shape[0],
            "n_samples": counts_df.shape[1] - 1,
        }

        del counts_df
        return {"meta": meta_entry, "gpl_results": gpl_results}


    max_workers = min(4, len(candidates))


    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:

        futures = [executor.submit(process_candidate, c) for c in candidates]

        for future in concurrent.futures.as_completed(futures):
            res = future.result()
            if res is not None:
                results_meta.append(res["meta"])
                merged_gpl_data.update(res["gpl_results"])

    st.session_state.gpl_data = merged_gpl_data

    gc.collect()

    return results_meta

def fetch_rnaseq_matrix(gse_id: str, meta):
    try:
        return fetch_and_normalize(
            accession=gse_id,
            gsm_ids=meta["gsm_ids"],
            geo_cache_dir="./geo_cache",
            save_output=False,
        )
    except Exception as e:
        st.error(f"geo_rnaseq_normalizer error: {e}")
        return None

# # # # # # # # # # # # # # # #
# MICROARRAY DATASET HANDLING #
# # # # # # # # # # # # # # # #

def _annotate_matrix(df: pd.DataFrame) -> pd.DataFrame:
    df = df.reset_index()
    df.rename(columns={df.columns[0]: "_probe_id"}, inplace=True)

    probe_ids = [str(pid).strip() for pid in df["_probe_id"]]
    looks_like_symbols = sum([int(is_gene_symbol(pid)) for pid in probe_ids[:min(100, len(probe_ids))]]) / 100

    if looks_like_symbols > 0.90:
        df.insert(0, "Name", df["_probe_id"])
        df = df.drop("_probe_id", axis=1)
        st.session_state.gpl_data = {}
        return df

    NULL_VALUES = {"---", "na", "", "null", "nan"}

    try:
        results = get_gene_symbol_column(
            df, st.session_state.selected_gpl, probe_ids=df["_probe_id"]
        )
    except ValueError as e:
        st.warning(f"{e} Keeping original probe IDs.")
        df.insert(0, "Name", df["_probe_id"])
        df = df.drop("_probe_id", axis=1)
        st.session_state.gpl_data = {}
        return df

    master_mapping: dict[str, str] = {}
    probe_series = df["_probe_id"].astype(str).str.strip()

    for gpl, item in results.items():
        gpl_table, matching_col, symbol_col = item[0], item[1], item[2]
        gpl_map = {}

        if symbol_col in gpl_table.columns and matching_col in gpl_table.columns:
            for probe, sym in zip(
                gpl_table[matching_col].astype(str).str.strip(),
                gpl_table[symbol_col].astype(str).str.strip(),
            ):
                if probe and sym and str(sym).lower() not in NULL_VALUES:
                    if "//" in sym:
                        sym = re.split(r'///|//', sym)[1].strip()
                    gpl_map[probe] = sym

        master_mapping.update(gpl_map)

    st.session_state.gpl_data = results
    df["Name"] = probe_series.map(master_mapping)

    before = len(df)
    df = df.dropna(subset=["Name"])
    after = len(df)
    if before != after:
        print(f"dropped {before - after} unmapped rows ({after} remaining)")
        log_record(
            "gene_annotation", "Probes with no gene symbol mapping were dropped.",
            dropped=before - after, remaining=after,
        )

    df = df.drop("_probe_id", axis=1)
    print(f"Remapped across {len(results)} GPL platform(s)")
    return df



def fetch_microarray_matrix(meta: dict):
    dp_text = st.session_state["dp_text"]
    norm_info = detect_normalization(dp_text=dp_text)

    log_record(
        "normalization", "Normalization detected from source metadata.",
        norm_type=norm_info.norm_type, is_log=norm_info.is_log, source=norm_info.source,
    )

    final_df = GLOBAL_GSE.pivot_samples(values="VALUE")[meta["gsm_ids"]]

    print(f"Microarray matrix contains {len(final_df.columns)} samples.")
    final_df = _annotate_matrix(final_df)
    final_df = reduce_matrix_memory(final_df)

    geo_cache_dir = Path("./geo_cache")
    geo_cache_dir.mkdir(parents=True, exist_ok=True)
    cache_id = uuid.uuid4().hex
    cache_path = geo_cache_dir / f"{gse_id}_{cache_id}_microarray.parquet"
    modified_path = geo_cache_dir / f"{gse_id}_{cache_id}_microarray_modified.parquet"

    # buffer = io.BytesIO()
    # final_df.to_parquet(buffer, index=False)
    # buffer.seek(0)

    # cache.add(cache_path, buffer.getvalue())
    # cache.add(modified_path, buffer.getvalue())

    final_df.to_parquet(cache_path)
    shutil.copy(cache_path, modified_path)

    meta_entry = {
        "path": cache_path,
        "modified_path": modified_path,
        "filename": "microarray_matrix",
        "normalization_type": norm_info.norm_type,
        "is_log": norm_info.is_log,
        "n_genes": final_df.shape[0],
        "n_samples": final_df.shape[1] - 1,
    }
    del final_df
    gc.collect()
    return [meta_entry]


def get_txt_gz_stream(parquet_path: str, columns):
    temp_dir = tempfile.gettempdir()
    output_path = os.path.join(temp_dir, f"export_{uuid.uuid4().hex}.txt.gz")

    if columns:
        if "Name" not in columns:
            columns = ["Name"] + columns
        cols_sql = ", ".join('"' + c.replace('"', '""') + '"' for c in columns)
        escaped_path = str(parquet_path).replace("'", "''")
        select_clause = f"SELECT {cols_sql} FROM '{escaped_path}'"
    else:
        escaped_path = str(parquet_path).replace("'", "''")
        select_clause = f"SELECT * FROM '{escaped_path}'"

    duckdb.execute(f"""
        COPY ({select_clause})
        TO '{output_path}'
        (FORMAT 'CSV', DELIMITER '\t', HEADER TRUE, COMPRESSION 'GZIP');
    """)

    try:
        with open(output_path, "rb") as f:
            data = f.read()
    finally:
        if os.path.exists(output_path):
            os.remove(output_path)

    return data


_NUMERIC_EXTRACT_RE = re.compile(r"[-+]?\d+(?:\.\d+)?")

def _extract_numeric(value) -> float:
    if pd.isna(value):
        return float("nan")
    s = str(value).strip()
    if not s:
        return float("nan")
    m = _NUMERIC_EXTRACT_RE.search(s)
    return float(m.group(0)) if m else float("nan")

_DEATH_WORDS = {"1", "1.0", "dead", "deceased", "died", "death", "event", "yes", "true"}
_ALIVE_WORDS = {"0", "0.0", "alive", "living", "censored", "no", "false", "no event"}

def _auto_binary_map(unique_vals) -> Optional[dict]:
    cleaned = {}
    for v in unique_vals:
        s = str(v).strip().lower()
        if s in _ALIVE_WORDS:
            cleaned[v] = 0
        elif s in _DEATH_WORDS:
            cleaned[v] = 1
        else:
            return None
    if len(unique_vals) == 2 and set(cleaned.values()) == {0, 1}:
        return cleaned
    return None

# # # # # # # # # # # # # # # #
# USER INTERFACE  # # # # # # #
# # # # # # # # # # # # # # # #

@st.fragment
def survival_metadata_ui(char_df):
    if char_df.empty or not len(char_df.columns):
        st.info("No sample characteristics are available for time-to-event data.")
        return

    if "survival_df" not in st.session_state:
        st.session_state.survival_df = pd.DataFrame(index=char_df.index)
        st.session_state.survival_df["GSM"] = char_df.index

    with st.expander("Generate Time-to-Event Data"):
        st.write("Preview of characteristics data (taken from NCBI GEO)")
        st.dataframe(
            char_df.head(),
            use_container_width=True,
            hide_index=True,
        )
        mortality = st.selectbox(
            "Event column",
            options=list(char_df.columns),
            key="survival_mortality_col"
        )

        # Iterate safely over unique values
        unique_vals = sorted(list(char_df[mortality].dropna().unique()))
        auto_map = _auto_binary_map(unique_vals)
        if auto_map is not None:
            st.session_state.survival_df["death"] = char_df[mortality].map(auto_map)
            st.caption(f"Auto-detected event coding: {auto_map}")
        else:
            mapping = {}
            for i in unique_vals:
                mapping[i] = st.radio(
                    label=f"Label '{i}' as no event (0) or event (1)",
                    options=(0, 1),
                    key=f"alive_dead_{i}"
                )
            st.session_state.survival_df["death"] = char_df[mortality].map(mapping)

        t_mortality = st.selectbox(
            "Time to event column",
            options=list(char_df.columns),
            key="survival_time_col"
        )

        time_units = st.radio(
            label=f"What are the units of time?",
            options=("days", "months", "years"),
        )

        st.session_state.survival_df.drop(
            columns=["days", "months", "years"],
            errors="ignore",
            inplace=True,
        )

        raw_time = char_df[t_mortality]
        parsed_time = raw_time.apply(_extract_numeric)
        n_unparseable = int(parsed_time.isna().sum() - raw_time.isna().sum())
        if n_unparseable > 0:
            st.warning(
                f"{n_unparseable} value(s) in '{t_mortality}' had no recognizable "
                f"number (e.g. 'unknown', 'N/A') and were set to missing."
            )
        st.session_state.survival_df[time_units] = parsed_time

        st.write("Final time-to-event data: ")

        st.dataframe(
            st.session_state.survival_df.head(),
            use_container_width=True,
            hide_index=True,
        )


        st.download_button(
            label="⬇ Download as .txt (tab-separated)",
            key = f"download_survival",
            data=st.session_state.survival_df.to_csv(sep="\t", index=False),
            file_name=f"{gse_id}_survival.txt",
            mime="text/plain",
        )

import os
import numpy as np
import pandas as pd
import plotly.express as px
import pyarrow.parquet as pq
import pyarrow.types as pt
import streamlit as st
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

MAX_ROWS_FOR_QC = 5000
SILHOUETTE_WARNING = 0.5


def _load_numeric(source_path):
    schema = get_schema(source_path)
    numeric_cols = [
        f.name for f in schema
        if pt.is_integer(f.type) or pt.is_floating(f.type)
    ]

    return get_parquet(source_path, columns = numeric_cols).to_pandas().astype(np.float32)


def _factor_labels(char_df, col, samples):
    lookup = char_df[col]
    lookup = lookup[~lookup.index.duplicated(keep="first")]
    return pd.Series(samples, index=samples).map(lookup).fillna("NA").astype(str)


@st.cache_data(show_spinner="PCA and RLE Loading")
def compute_microarray_qc(source_path, char_df, raw = False):
    source_df = _load_numeric(source_path)
    if raw:
            source_df = source_df[(source_df >= 10).sum(axis=1) >= 3]
            cpm = source_df / source_df.sum(axis=0) * 1e6
            logcpm = np.log2(cpm + 1)
            top = logcpm.var(axis=1).nlargest(1000).index
            X = logcpm.loc[top].T

    if source_df.max().max() > 100:
        source_df = np.log2(source_df.clip(lower=0) + 1)

        st.warning("Existing data was likely linearly scaled, so for PCA and RLE purposes it has been log scaled.")

    if not raw:
        X = source_df.T
        X = X.fillna(X.mean()).dropna(axis=1)

    pca = PCA(n_components=2, svd_solver="randomized", random_state=0)
    pcs = pca.fit_transform(X)
    pca_base = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=X.index)
    pca_base["Sample"] = pca_base.index
    var_explained = pca.explained_variance_ratio_ * 100

    if len(source_df) > MAX_ROWS_FOR_QC:
        source_df = source_df.sample(MAX_ROWS_FOR_QC, random_state=0)
    deviation_df = source_df.sub(source_df.median(axis=1), axis=0)
    box_df = deviation_df.melt(var_name="Sample", value_name="Deviation")

    silhouettes = {}
    for col in char_df.columns:
        labels = _factor_labels(char_df, col, pca_base["Sample"])
        if 1 < labels.nunique() < len(labels):
            silhouettes[col] = silhouette_score(pcs, labels)

    return pca_base, var_explained, box_df, silhouettes


@st.cache_data(show_spinner="Sequencing Depth Loading ")
def compute_depth(source_path):
    source_df = _load_numeric(source_path)
    return pd.DataFrame({"Sample": source_df.columns, "Depth": source_df.sum(axis=0).values})




def qc(char_df, gse_id, source_path, modified_path, normalization_type, selected_columns):
    with st.expander("Quality Control", expanded=False):
        st.session_state["selected_columns_for_download"] = selected_columns
        selected = set(selected_columns)

        depth_df = None

        if "microarray" in str(source_path):
            pca_base, var_explained, box_df, silhouettes = compute_microarray_qc(
                str(source_path), char_df
            )

            best_col = max(silhouettes, key=silhouettes.get) if silhouettes else char_df.columns[0]

            st.markdown(
                "#### Sample PCA\n"
                "Each point is one sample. Samples with similar gene expression sit close together. "
                "Use this to spot **outliers** (points far from everything else) that may be bad "
                "samples worth removing."
            )
            color_col = st.selectbox(
                "Varied Factor (color the samples by...)",
                char_df.columns,
                list(char_df.columns).index(best_col),
                key=f"qc_rle_selector",
                help="Defaults to the factor that best separates samples on the PCA.",
            )

            pca_df = pca_base.copy()
            pca_df[color_col] = _factor_labels(char_df, color_col, pca_df["Sample"]).values
            pca_df["Status"] = np.where(
                pca_df["Sample"].isin(selected), "Kept", "Removed"
            )

            score = silhouettes.get(color_col)
            if score is not None and score > SILHOUETTE_WARNING:
                st.warning(
                    f"Samples separate very strongly by **{color_col}** "
                    f"(silhouette score = {score:.2f}). If this factor is not a real biological "
                    "difference (e.g. a processing date, platform, or lab), there is a strong "
                    "potential for **batch effects**. Consider correcting for it before "
                    "comparing groups."
                )

            n_removed = int((pca_df["Status"] == "Removed").sum())
            st.caption(
                f"Circle = sample kept (will be in your download)   |   "
                f"X = sample removed (excluded from your download)   |   "
                f"{len(pca_df) - n_removed} kept, {n_removed} removed"
                + (f"   |   Silhouette for '{color_col}': {score:.2f} "
                   "(-1 to 1; higher = groups are more distinct)" if score is not None else "")
            )


            pca_fig = px.scatter(
                pca_df,
                x="PC1",
                y="PC2",
                color=color_col,
                symbol="Status",
                symbol_map={"Kept": "circle", "Removed": "x"},
                hover_name="Sample",
                labels={
                    "PC1": f"PC1 ({var_explained[0]:.1f}% of variance)",
                    "PC2": f"PC2 ({var_explained[1]:.1f}% of variance)",
                },
                title="Sample PCA",
            )
            pca_fig.update_traces(marker=dict(size=10))
            st.plotly_chart(pca_fig, key=f"qc_pca")

            st.markdown(
                "#### Relative Log Expression (RLE)\n"
                "Each box shows how far one sample's genes deviate from the typical (median) "
                "gene value of each gene in the sample. In a healthy dataset, boxes are **centered near 0 and roughly the "
                "same height**. A box that is shifted up/down or much wider than the rest may "
                "be a low-quality sample or a batch difference."
            )
            box_df = box_df.copy()
            box_df[color_col] = box_df["Sample"].map(
                dict(zip(pca_df["Sample"], pca_df[color_col]))
            )
            box_df = box_df.sort_values(by=color_col)

            fig = px.box(box_df, x="Sample", y="Deviation", color=color_col, points=False)
            st.plotly_chart(fig, key=f"qc_plot")

            qc_summary_df = pca_df[["Sample", color_col, "Status"]].rename(
                columns={color_col: "colored_by"}
            )
            if normalization_type == "raw_counts" and depth_df is not None:
                qc_summary_df = qc_summary_df.merge(
                    depth_df[["Sample", "Depth"]], on="Sample", how="left"
                )

            st.download_button(
                label="⬇ Download QC summary (.csv)",
                data=qc_summary_df.to_csv(index=False),
                file_name=f"{gse_id}_qc_summary.csv",
                mime="text/csv",
                key="download_qc_summary",
            )

        else:
            if normalization_type == "raw_counts":
                pca_base, var_explained, box_df, silhouettes = compute_microarray_qc(
                    str(source_path), char_df, raw=True
                )
            else:
                pca_base, var_explained, box_df, silhouettes = compute_microarray_qc(
                    str(source_path), char_df, raw=True
                )

            best_col = max(silhouettes, key=silhouettes.get) if silhouettes else char_df.columns[0]

            st.markdown(
                "#### Sample PCA\n"
                "Each point is one sample. Samples with similar gene expression sit close together. "
                "Use this to spot **outliers** (points far from everything else) that may be bad "
                "samples worth removing."
            )
            color_col = st.selectbox(
                "Varied Factor (color the samples by...)",
                char_df.columns,
                list(char_df.columns).index(best_col),
                key=f"qc_rle_selector",
                help="Defaults to the factor that best separates samples on the PCA.",
            )

            pca_df = pca_base.copy()
            pca_df[color_col] = _factor_labels(char_df, color_col, pca_df["Sample"]).values
            pca_df["Status"] = np.where(
                pca_df["Sample"].isin(selected), "Kept", "Removed"
            )

            score = silhouettes.get(color_col)
            if score is not None and score > SILHOUETTE_WARNING:
                st.warning(
                    f"Samples separate very strongly by **{color_col}** "
                    f"(silhouette score = {score:.2f}). If this factor is not a real biological "
                    "difference (e.g. a processing date, platform, or lab), there is a strong "
                    "potential for **batch effects**. Consider correcting for it before "
                    "comparing groups."
                )

            n_removed = int((pca_df["Status"] == "Removed").sum())
            st.caption(
                f"Circle = sample kept (will be in your download)   |   "
                f"X = sample removed (excluded from your download)   |   "
                f"{len(pca_df) - n_removed} kept, {n_removed} removed"
                + (f"   |   Silhouette for '{color_col}': {score:.2f} "
                   "(-1 to 1; higher = groups are more distinct)" if score is not None else "")
            )


            pca_fig = px.scatter(
                pca_df,
                x="PC1",
                y="PC2",
                color=color_col,
                symbol="Status",
                symbol_map={"Kept": "circle", "Removed": "x"},
                hover_name="Sample",
                labels={
                    "PC1": f"PC1 ({var_explained[0]:.1f}% of variance)",
                    "PC2": f"PC2 ({var_explained[1]:.1f}% of variance)",
                },
                title="Sample PCA",
            )
            pca_fig.update_traces(marker=dict(size=10))
            st.plotly_chart(pca_fig, key=f"qc_pca")



            if normalization_type == "raw_counts":
                depth_df = compute_depth(str(source_path))
                depth_df["Status"] = np.where(
                    depth_df["Sample"].isin(selected), "Kept", "Removed"
                )
                st.markdown(
                    "#### Sequencing Depth\n"
                    "Total reads per sample. Samples below the dashed line (20 million reads) have "
                    "less data, so their gene measurements are noisier and less reliable. "
                    "Works only for raw counts."
                )
                fig = px.bar(
                    depth_df, x="Sample", y="Depth",
                    labels={"Depth": "Total reads"},
                    color = "Status",
                    color_discrete_map = {
                            "Kept": "blue",
                            "Removed": "yellow"
                        }

                )
                fig.add_hline(
                    y=20e6,
                    line_dash="dash",
                    line_color="black",
                    annotation_text="20M reads",
                )
                st.plotly_chart(fig, key=f"qc_plot")




@st.fragment
def annotate_columns(i, char_df, gse_id, source_path, modified_path):
    with st.expander("Annotate Columns", expanded=False):
        st.write("Replace GSM column names with sample data")

        st.pills("Annotate by:", char_df.columns, selection_mode="multi", key=f"pill_selector_{i}")

        current_selection = st.session_state.get(f"pill_selector_{i}") or []
        if current_selection:
            example_name = "_".join([str(char_df[x].iloc[0]) for x in current_selection])
            st.write(f"Example sample name: {example_name}_{char_df.index[0]}")
        else:
            st.write(char_df.index[0])

        if st.button("Annotate Columns", key=f"annotate_columns_{i}"):
            column_mapping = {}
            if not current_selection:
                for gsm in char_df.index:
                    column_mapping[gsm] = str(gsm)
            else:
                for gsm in char_df.index:
                    column_mapping[gsm] = "_".join([str(char_df[x].loc[gsm]) for x in current_selection]) + f"_{gsm}"

            modified_df = get_parquet(modified_path)

            dynamic_mapping = {}
            for col in modified_df.column_names:
                for k, v in column_mapping.items():
                    if k in col:
                        dynamic_mapping[col] = v
                        break

            new_column_names = [
                dynamic_mapping.get(col, col) for col in modified_df.column_names
            ]

            modified_df = modified_df.rename_columns(new_column_names)

            log_record(
                "sample_annotation", "GSM column names replaced with characteristic-derived names.",
                fields_used=current_selection or "GSM ID only",
            )


            pq.write_table(
                modified_df,
                modified_path,
                compression="snappy"
            )

            st.session_state[f"selected_cols_dict_{gse_id}_{i}"] = {col: True for col in list(column_mapping.values())}
            st.session_state[f"selected_cols_dict_{gse_id}_{i}"]["Name"] = True
            st.rerun(scope="app")

@st.fragment
def column_selector(counts_path, i, gse_id, state, suffix):
    parquet_file = pq.ParquetFile(counts_path)
    all_columns = parquet_file.schema.names

    state_key = f"selected_cols_dict_{gse_id}_{state}_{suffix}"
    if state_key not in st.session_state:
        st.session_state[state_key] = {col: True for col in all_columns}

    st.markdown("### Select Columns to Include in Download \n" \
    "You can search for specific columns using the search bar below. Uncheck any columns you wish to exclude from your download." \
    "You can use the Quality Control data below to inform your column selections.")
    search_term = st.text_input("Search columns", key=f"search_bar_{i}_{suffix}")

    filtered_columns = [
        col for col in all_columns
        if search_term.lower() in col.lower()
    ] if search_term else list(all_columns)

    def select_all():
        for col in filtered_columns:
            st.session_state[state_key][col] = True
            st.session_state[f"ui_{col}_{i}_{suffix}_key"] = True
        st.session_state["_force_full_rerun"] = True

    def deselect_all():
        for col in filtered_columns:
            st.session_state[state_key][col] = False
            st.session_state[f"ui_{col}_{i}_{suffix}_key"] = False
        st.session_state["_force_full_rerun"] = True

    col1, col2 = st.columns(2)
    with col1:
        st.button("Select All", on_click=select_all, key=f"select_all_{i}_{suffix}")
    with col2:
        st.button("Deselect All", on_click=deselect_all, key=f"deselect_all_{i}_{suffix}")

    with st.expander("Show All Columns", expanded=False):
        num_grid_cols = 4
        grid_columns = st.columns(num_grid_cols)

        for idx, col in enumerate(filtered_columns):
            def update_single_col(c=col):
                st.session_state[state_key][c] = st.session_state[f"ui_{c}_{i}_{suffix}_key"]
                st.session_state["_force_full_rerun"] = True


            widget_key = f"ui_{col}_{i}_{suffix}_key"

            if col not in st.session_state[state_key]:
                st.session_state[state_key][col] = True

            st.session_state[widget_key] = st.session_state[state_key][col]

            with grid_columns[idx % num_grid_cols]:
                st.checkbox(
                    col,
                    key=widget_key,
                    on_change=update_single_col
                )

    selected_columns = [
        col for col in all_columns
        if st.session_state[state_key].get(col, True)
    ]

    return selected_columns

def get_parquet(path, **kwargs):
    # reader = pa.BufferReader(cache.get(path))
    # return pq.read_table(reader, **kwargs)
    return pq.read_table(path, **kwargs)

def get_schema(path, **kwargs):
    # reader = pa.BufferReader(cache.get(path))
    # return pq.read_schema(reader, **kwargs)
    return pq.read_schema(path, **kwargs)

def apply_norm(modified_path, norm, col_list, i, source_path):
    target_cols = set(col_list)

    schema = get_schema(modified_path)
    base_schema = get_schema(source_path)
    base_column_set = set(base_schema.names)

    numeric_cols = [
        field.name for field in schema
        if (pt.is_integer(field.type) or pt.is_floating(field.type)) and field.name in target_cols
    ]

    base_gsm_cols = []
    valid_numeric_cols = []

    for col in numeric_cols:
        gsm_id = col.rsplit("_", 1)[-1]
        if gsm_id in base_column_set:
            base_gsm_cols.append(gsm_id)
            valid_numeric_cols.append(col)
        elif col in base_column_set:
            base_gsm_cols.append(col)
            valid_numeric_cols.append(col)

    if not base_gsm_cols:
        return get_parquet(modified_path)

    numeric_df = get_parquet(source_path, columns = base_gsm_cols)
    numeric_df = numeric_df.rename_columns(valid_numeric_cols)

    if norm == "none":
        new_cols_df = numeric_df

    elif norm == "log2":
        new_cols_df = pa.Table.from_arrays(
            [pc.log2(pc.add(numeric_df.column(col), 1.0)) for col in valid_numeric_cols],
            names=valid_numeric_cols
        )

    elif norm == "log10":
        new_cols_df = pa.Table.from_arrays(
            [pc.log10(pc.add(numeric_df.column(col), 1.0)) for col in valid_numeric_cols],
            names=valid_numeric_cols
        )

    elif norm == "cpm":
        processed_arrays = []
        for col in valid_numeric_cols:
            col_arr = numeric_df.column(col)
            col_sum = pc.sum(col_arr).as_py()

            cpm_arr = pc.multiply(pc.divide(col_arr, col_sum), 1e6)
            processed_arrays.append(cpm_arr)

        new_cols_df = pa.Table.from_arrays(processed_arrays, names=valid_numeric_cols)

    elif norm == "log2(cpm+1)":
        processed_arrays = []
        for col in valid_numeric_cols:
            col_arr = numeric_df.column(col)
            col_sum = pc.sum(col_arr).as_py()
            cpm_arr = pc.multiply(pc.divide(col_arr, col_sum), 1e6)

            log_cpm = pc.log2(pc.add(cpm_arr, 1.0))
            processed_arrays.append(log_cpm)

        new_cols_df = pa.Table.from_arrays(processed_arrays, names=valid_numeric_cols)

    elif norm == "log10(cpm+1)":
        processed_arrays = []
        for col in valid_numeric_cols:
            col_arr = numeric_df.column(col)
            col_sum = pc.sum(col_arr).as_py()
            cpm_arr = pc.multiply(pc.divide(col_arr, col_sum), 1e6)

            log_cpm = pc.log10(pc.add(cpm_arr, 1.0))
            processed_arrays.append(log_cpm)

        new_cols_df = pa.Table.from_arrays(processed_arrays, names=valid_numeric_cols)

    else:
        del numeric_df
        gc.collect()
        return get_parquet(modified_path)

    del numeric_df
    gc.collect()

    df = get_parquet(modified_path)

    overlapping_cols = set(valid_numeric_cols).intersection(df.column_names)
    if overlapping_cols:
        df = df.drop_columns(list(overlapping_cols))

    for col_name in valid_numeric_cols:
        df = df.append_column(col_name, new_cols_df.column(col_name))

    del new_cols_df
    gc.collect()

    return df

if "result_lists" not in st.session_state:
    st.session_state.result_lists = None

@st.dialog("Download Data", width="large")
def download_dialog():
    if st.session_state.result_lists is not None:
        result_lists = st.session_state.result_lists
        current_norms = []
        for i in range(len(result_lists)):
            current_norms.append(st.session_state.get(f"norm_type_{i}", "none"))
        def prep_download(i, meta_entry, columns = None):
            current_norm = current_norms[i]
            norm_suffix = f"_{current_norm}" if current_norm != "none" else ""
            return (f"download_{i}", get_txt_gz_stream((meta_entry["modified_path"]),
                                                        columns = columns), f"{gse_id}_{meta_entry['normalization_type']}{norm_suffix}.txt.gz",
                                                        meta_entry)

        select_dict = {}
        for i, meta_entry in enumerate(result_lists):
            select_dict[meta_entry["modified_path"]] = (i, meta_entry)

        selected_path = st.selectbox(
            "Which dataset would you like to download?",
            options = [x["modified_path"] for x in result_lists],
            key = "selected_path_for_download"
        )

        selected_columns = column_selector(
            selected_path, "download", gse_id, "qc", selected_path
        )

        char_df = st.session_state.char_df.copy()
        gsm_filter = []
        for g in list(st.session_state.current_gpl.keys()):
            if g != "ncbi":
                gsm_filter.extend(meta["gsm_gpl_dict"][g])
        if not gsm_filter:
            gsm_filter = meta["gsm_ids"]


        res = prep_download(select_dict[st.session_state["selected_path_for_download"]][0], select_dict[st.session_state["selected_path_for_download"]][1],
                            columns = selected_columns)
        if res is not None:
            # signature = st.text_area(
            #     "Paste signature to check availability (any punctuation okay)"
            # )

            # signature_list = re.split(r'[;,\s]', signature)

            # signature_set = set([s.strip() for s in signature_list if s.strip()])


            # if signature_set:
            #     hngc_data = fetch_hgnc_table()
            #     df = pq.read_table(meta_entry["modified_path"], columns=["Name"]).to_pandas()
            #     available_genes = set(df["Name"].dropna().unique())
            #     print(list(available_genes)[:10])

            #     signature_set_symbol = signature_set.intersection(set(hngc_data["symbol"].dropna().unique()))
            #     signature_set_alias = signature_set.intersection(set(hngc_data["alias_symbol"].dropna().unique()))
            #     signature_set_prev = signature_set.intersection(set(hngc_data["prev_symbol"].dropna().unique()))
            #     missing_genes_symbol = signature_set_symbol - available_genes
            #     missing_genes_alias = signature_set_alias - available_genes
            #     missing_genes_prev = signature_set_prev - available_genes

            #     print(missing_genes_symbol.intersection(missing_genes_alias).intersection(missing_genes_prev))

            #     if intersection := missing_genes_symbol.intersection(missing_genes_alias).intersection(missing_genes_prev):

            #         st.warning(
            #             f"{len(intersection)} gene(s) from your signature were not found in this dataset: "
            #             + ", ".join(list(intersection)[:10]) + ("..." if len(intersection) > 10 else "")
            #         )

            #         log_record(
            #             "signature_check", "User checked signature against dataset.",
            #             missing_genes=list(intersection),
            #         )
            #     else:
            #         st.success("All genes from your signature are present in this dataset.")

            qc_meta_entry = res[3]
            qc(char_df.loc[char_df.index.isin(gsm_filter) | ~char_df.index.str.startswith("GSM")], gse_id, qc_meta_entry["path"], qc_meta_entry["modified_path"], qc_meta_entry["normalization_type"], selected_columns)
            if st.session_state.get("_force_full_rerun"):
                st.session_state["_force_full_rerun"] = False
                st.rerun()

            st.markdown("### Download Dataset \n" \
            "Click the button below to download the dataset as a tab-separated text file. The file will be compressed using gzip to reduce download size.")

            st.download_button(
                label="Download Dataset ⬇",
                key=res[0],
                data=res[1],
                file_name=res[2],
                mime="application/gzip",
            )

        if st.session_state.get("record_log"):

            st.markdown("### Download Provenance Record \n" \
            "Click the button below to download a JSON file containing a record of the actions taken during the analysis. This can be used for reproducibility and tracking of the analysis steps.")

            record_payload = {
                "gse_id": gse_id,
                "accessed": datetime.now(timezone.utc).isoformat(),
                "tool_version": "geo-2-compass-0.1",
                "events": st.session_state.record_log,
            }
            st.download_button(
                label="⬇ Download provenance (.json)",
                data=json.dumps(record_payload, indent=2, default=str),
                file_name=f"{gse_id}_record.json",
                mime="application/json",
                key="download_record",
            )


    else:
        st.info("No results available yet.")

c1, c2 = st.columns(2)
with c1:
    if st.button("Fetch & Build Matrix", type="primary"):
        st.session_state.run_pipeline = True
        st.session_state.result_lists = None

        for key in list(st.session_state.keys()):
            if str(key).startswith(("df_", "base_df_", "norm_type_", "preview_", "selected_cols_dict_", "pill_selector_")):
                del st.session_state[key]
        gc.collect()
with c2:
    if st.button("Download Data", type = "primary"):
        download_dialog()



if st.session_state.run_pipeline and st.session_state.result_lists is None:

    with st.status("Running pipeline…", expanded=True) as status:

        if meta["type"] == "Microarray":
            st.write("Running microarray pipeline (async GSM fetch + GPL annotation)…")
            result_lists = fetch_microarray_matrix(meta)
        else:
            st.write(
                f"Running RNA-seq pipeline for **{gse_id}** "
            )
            result_lists = fetch_rnaseq_matrix(gse_id, meta)

        if not result_lists:
            status.update(label="Pipeline failed.", state="error")
            st.error("No usable expression matrix was found for this accession.")
            st.stop()
        st.session_state.result_lists = result_lists
        status.update(label="Done!", state="complete")

if st.session_state.result_lists is not None:
    result_lists = st.session_state.result_lists

    st.session_state.current_gpl = select_gpl(st.session_state.gpl_data)

    char_df = st.session_state.char_df.copy()


    gsm_filter = []

    for g in list(st.session_state.current_gpl.keys()):
        if g != "ncbi":
            gsm_filter.extend(meta["gsm_gpl_dict"][g])

    if not gsm_filter:
        gsm_filter = meta["gsm_ids"]

    survival_metadata_ui(char_df.loc[char_df.index.isin(gsm_filter) | ~char_df.index.str.startswith("GSM")])


    qc_meta_entry = result_lists[0]

    for i, meta_entry in enumerate(result_lists):
        if meta_entry["normalization_type"] == "raw_counts":
            qc_meta_entry = meta_entry

    # qc(char_df.loc[char_df.index.isin(gsm_filter) | ~char_df.index.str.startswith("GSM")], gse_id, qc_meta_entry["path"], qc_meta_entry["modified_path"], qc_meta_entry["normalization_type"])
    st.markdown("""
        <style>
            .normalization {
                font-size: 16px !important;
                background-color: #FFCCCC; /* Light red */
                border: 2px rgb(255, 75, 75) solid;
                text-align: center;
                color: #000000;
                padding: 0.25rem 0.75rem;
                margin-left: 10%;
                margin-right: 10%;
                border-radius: 0.375rem 0.375rem 0 0;
            }

            .container {
                padding: 0;
                border-radius: 0.5rem;
                border: 2px rgba(100, 100, 100, 0.2) solid;
            }

            .normtext {
                padding: 0.5rem;
                font-size: 16px !important;
                text-align: center;
            }

        </style>
        """, unsafe_allow_html=True)

    for i, meta_entry in enumerate(result_lists):
        if f"norm_type_{i}" not in st.session_state:
            st.session_state[f"norm_type_{i}"] = "none"


        original_normalization = meta_entry["normalization_type"]
        n_genes = meta_entry["n_genes"]
        n_samples = meta_entry["n_samples"]

        st.markdown(body="<hr>", unsafe_allow_html=True)
        st.caption(f"**{meta_entry['filename']}**")
        c4, c5 = st.columns(2)
        with c4:
            st.markdown(body=f'<div class="container"><p class="normalization">Original Normalization</p> <p class="normtext">{original_normalization}</p></div>', unsafe_allow_html=True)
        with c5:
            st.markdown(body=f'<div class="container"><p class="normalization">Applied Normalization</p> <p class="normtext">{st.session_state[f"norm_type_{i}"].upper()}</p></div>', unsafe_allow_html=True)

        c1, c2, c3 = st.columns(3)
        c1.metric("Genes", f"{n_genes:,}")
        c2.metric("Samples", len(char_df.loc[char_df.index.isin(gsm_filter) | ~char_df.index.str.startswith("GSM")]))

        print(gsm_filter)

        annotate_columns(i, char_df.loc[char_df.index.isin(gsm_filter) | ~char_df.index.str.startswith("GSM")], gse_id, meta_entry["path"], meta_entry["modified_path"])

        with st.expander("Change Normalization", expanded=False):
            st.write(st.session_state["dp_text"])
            is_log = meta_entry.get("is_log", False)
            st.write("What normalization would you like to apply?")
            if original_normalization == "raw_counts":
                norm_suggestion = f"log2 or log2(cpm+1) because {original_normalization} data is raw."
            elif is_log:
                norm_suggestion = f"None because {original_normalization} already log-scales the data."
            else:
                norm_suggestion = f"log2 because {original_normalization} is linearly scaled, potentially leading to a skewed distribution."

            st.radio(
                label=f"Suggestion: {norm_suggestion}",
                options=("none", "log2", "log10", "cpm", "log2(cpm+1)", "log10(cpm+1)"),
                key=f"norm_type_{i}"
            )

            st.write("Which columns would you like to apply it to?")
            st.session_state[f"selected_columns_{i}"] = column_selector(meta_entry["modified_path"], i, gse_id, "annotator", str(meta_entry["modified_path"]) + "&")

            submit_button = st.button(label="Renormalize", key=f"renormalize_{i}")

            if submit_button:
                modified_df = apply_norm(
                    meta_entry["modified_path"],
                    st.session_state[f"norm_type_{i}"],
                    st.session_state[f"selected_columns_{i}"],
                    i,
                    meta_entry["path"]
                )

                log_record(
                    "normalization", "User applied a renormalization.",
                    applied=st.session_state[f"norm_type_{i}"],
                    original=original_normalization,
                    n_columns=len(st.session_state[f"selected_columns_{i}"]),
                )

                pq.write_table(
                    modified_df,
                    meta_entry["modified_path"],
                    compression="snappy"
                )



        # Live view of the CURRENT df_key — reflects annotate/renormalize immediately
        parquet_file = pq.ParquetFile(meta_entry["modified_path"])
        total_rows = parquet_file.metadata.num_rows
        schema_names = parquet_file.schema.names

        gsm_filter = {
            gsm
            for g in st.session_state.current_gpl.keys()
            for gsm in meta.get("gsm_gpl_dict", {}).get(g, [])
        }

        if gsm_filter == set():
            gsm_filter = set(meta.get("gsm_ids", []))

        all_cols = []
        for col in schema_names:
            if "GSM" not in col.upper():
                all_cols.append(col)
            else:
                match = re.search(r"GSM\d+", col, re.IGNORECASE)
                if match and match.group(0).upper() in gsm_filter:
                    all_cols.append(col)
        total_cols = len(all_cols)

        if "Name" in all_cols[:PREVIEW_COLUMNS]:
            preview_cols = all_cols[:PREVIEW_COLUMNS]
        else:
            preview_cols = ["Name"] + all_cols[:PREVIEW_COLUMNS-1]

        if total_rows > PREVIEW_ROWS or total_cols > PREVIEW_COLUMNS:
            st.caption(
                f"Showing {min(total_rows, PREVIEW_ROWS):,}/{total_rows:,} rows and "
                f"{len(preview_cols):,}/{total_cols:,} columns. The download contains the complete matrix."
            )

        dataset = ds.dataset(meta_entry["modified_path"], format="parquet")
        preview_table = dataset.head(PREVIEW_ROWS, columns=preview_cols)

        st.dataframe(
            preview_table.to_pandas(),
            use_container_width=True,
            hide_index=True,
            column_config={"Name": st.column_config.TextColumn("Name", pinned=True, width="small")},
        )


        current_norm = st.session_state.get(f"norm_type_{i}", "none")
        norm_suffix = f"_{current_norm}" if current_norm != "none" else ""

        # st.download_button(
        #     label="⬇ Download complete matrix (.txt.gz)",
        #     key=f"download_{i}",
        #     data=StreamWrapper(stream_parquet_to_tsv_gz(meta_entry["modified_path"])),
        #     file_name=f"{gse_id}_{original_normalization}{norm_suffix}.txt.gz",
        #     mime="application/gzip",
        # )
