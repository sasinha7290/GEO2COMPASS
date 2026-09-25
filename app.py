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

# PyArrow 25.0.0 can segfault when Streamlit initializes Arrow from a
# ScriptRunner thread. Use the system allocator even if the deployment
# environment does not define this variable .
os.environ.setdefault("ARROW_DEFAULT_MEMORY_POOL", "system")

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.types as pt
import pyarrow.dataset as ds
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

MAX_DOWNLOAD_BYTES = 200 * 1024 * 1024
MAX_MATRIX_CELLS = 50_000_000
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
        old_gse_id = st.session_state.get("active_gse_id")
        if old_gse_id:
            geo_cache_dir = Path("./geo_cache")
            for f in geo_cache_dir.glob(f"{old_gse_id}_*.parquet"):
                f.unlink(missing_ok=True)
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
                delete=True,
            ) as temp_file:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        temp_file.write(chunk)

                temp_file.flush()
                print("Parsing data with GEOparse from temporary file...")
                return GEOparse.get_GEO(
                    filepath=temp_file.name,
                    geotype="GSE",
                )
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Download failed: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Failed to parse GEO data: {e}") from e

def get_dp_and_char(gse, meta):
    first_gsm = gse.gsms[meta["gsm_ids"][0]]
    dp_text = first_gsm.metadata.get("data_processing", [""])[0]
    st.session_state["dp_text"] = dp_text
    char_list = list(dict.fromkeys(
        x.split(": ", 1)[0]
        for x in first_gsm.metadata.get("characteristics_ch1", [])
    ))
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

def needs_log(dp_text: str):
    yes_log = ["rma ", "(rma)", "lowess", "log", "vsn", "beadstudio","vsn","fhma","plier","quantile"]
    no_log = ["mas5", "mas 5"]
    for n in no_log:
        if n in dp_text.lower():
            return False, n
    for y in yes_log:
        if y in dp_text.lower(): 
            return True, y
    
    if "raw" in dp_text.lower():
        return False, "raw"
    
    return False, dp_text

def to_cpm(df):
    numeric_df = df.apply(pd.to_numeric, downcast = "float", errors="coerce")
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
        gpl_df["new_id_col"] = gpl_df[best_col].astype(str).str.replace(r"_at$", "", regex=True)
        if (str(gpl_df[best_col].iloc[0]).strip()[:3].lower() == "eg:"):
            gpl_df["new_id_col"] = gpl_df[best_col].astype(str).str.split(':').str[1]

        try:
            gp = GProfiler(return_dataframe=True)
            results = gp.convert(organism=species, query=list(gpl_df["new_id_col"]),
                                numeric_namespace=numeric_namespace, target_namespace=target_namespace)
            if results.empty or "name" not in results.columns:
                #st.warning("Gene symbol conversion returned no results; keeping original IDs.")
                return gpl_df, symbol_col
        except Exception as e:
            #st.warning(f"Gene symbol conversion failed ({e}); keeping original IDs.")
            return gpl_df, symbol_col
        
        results_deduped = results.drop_duplicates(subset="incoming", keep="first")
        mapping = results_deduped.set_index("incoming")["name"]
        gpl_df["gene_symbol"] = gpl_df["new_id_col"].values 
        gpl_df["gene_symbol"] = [mapping.get(id_, float("nan")) for id_ in gpl_df["new_id_col"]]

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
            print("LOLLAPALOOZA")

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
    

#normlization

def classify_normalization(dp_text: str) -> str:
    norm_type = "unknown"
    if (any(s in dp_text.lower() for s in ["raw","rsem", "unnorm","htseq-count","feature"])):
        norm_type = "raw_counts"
    elif (any(s in dp_text.lower() for s in ["cpm"])):
        norm_type = "cpm"
    elif (any(s in dp_text.lower() for s in ["tpm"])):
        norm_type = "tpm"
    elif (any(s in dp_text.lower() for s in ["rpkm","fpkm"])):
        norm_type = "rpkm_fpkm"
    elif (any(s in dp_text.lower() for s in ["quantile"])):
        norm_type = "quantile"
    elif (any(s in dp_text.lower() for s in ["geometric"])):
        norm_type = "geometric"
    else:
        norm_type = "unknown"
        print("No clear normalization found.")
    return norm_type
 
def _download_bytes(url: str, cap = MAX_DOWNLOAD_BYTES) -> Optional[bytes]:
    buf = io.BytesIO()
    with requests.get(url, stream=True, timeout=(30, 300)) as r:
        r.raise_for_status()
        for chunk in r.iter_content(1024 * 1024):
            buf.write(chunk)
            if buf.tell() > cap:
                return None
    return buf.getvalue()

 
 
def _download_all(urls: list[str]) -> dict[str, Optional[bytes]]:
    results: dict[str, Optional[bytes]] = {}
    if not urls:
        return results
 
    # Limit workers so we don't spam open thousands of connections simultaneously
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
                    print(f"    ✓ {fname} ({len(res) / 1e6:.1f} MB)")
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

        if (any(s in lower for s in ["raw", "rsem", "unnorm","htseq-count","feature"])):
            normalization = "raw_counts"
        elif "tpm" in lower:
            normalization = "tpm"
        elif (any(s in lower for s in ["rpkm","fpkm"])):
            normalization = "fpkm_rpkm"
        elif "cpm" in lower:
            normalization = "cpm"
        elif (any(s in lower for s in ["quantile"])):
            normalization = "quantile"
        elif (any(s in lower for s in ["geometric"])):
            normalization = "geometric"
        else:
            normalization = "unknown"

        print(f"\n\n\n {filename} \t {normalization} \n\n\n")
            
 
        if not (is_tar or _is_tabular(filename)) or "annot" in filename:
            print(f"    [skip] {filename}")
            continue
        
        files.append(FileMeta(url=full_url, filename=filename, is_tar=is_tar, ncbi_data=ncbi_data, normalization=normalization))
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
        raw = url_bytes.get(meta.url)
        if raw is None:
            continue
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tf:
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
        raw = url_bytes.get(meta.url)
        if raw is None:
            continue
        filename = meta.filename
        try:
            if filename.endswith(".gz"):
                raw = gzip.decompress(raw)
                filename = filename[:-3]

            lower = filename.lower()
            if lower.endswith(".xls") or lower.endswith(".xlsx"):
                df = pd.read_excel(io.BytesIO(raw), index_col=0)
            else:
                sr = detect_delimiter(raw)
                df = pd.read_csv(
                    io.BytesIO(raw), sep=sr[0],
                    lineterminator= sr[1],
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
    inferred_norm = classify_normalization(dp_text)
    if not candidates:
        return []

    results_meta = []
    merged_gpl_data = {}

    # url_bytes = _download_all([c.url for c in candidates])

    #annot_url = "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
    #annot_bytes = _download_all([annot_url])[annot_url]
    #raw = url_bytes.get(annot_url)

    all_selected_gpls = st.session_state.selected_gpl

    def process_candidate(c):
        selected = [c]
        if c.normalization == "unknown":
            c.normalization = inferred_norm

        print(f"\n  Processing: {c.filename}")

        url_bytes = _download_all([c.url])

        if url_bytes.get(c.url) is None:
            print(f"    Exceeded file size, skipping {c.url}")
            return None
        
        counts_df = _fetch_counts_df(accession, selected, url_bytes)

        gc.collect()

        if counts_df.empty:
            return None

        if counts_df.size > MAX_MATRIX_CELLS:
            print(f"    Skipping {c.filename}: {counts_df.size:,} cells exceeds "
                  f"{MAX_MATRIX_CELLS:,}-cell limit")
            return None

        counts_df, gpl_results = _annotate_counts(accession, selected, counts_df, all_selected_gpls)
        #counts_df = reduce_matrix_memory(counts_df)

        retained_cols = [
            col for col in counts_df.columns
            if "gsm" not in col.lower() or any(gsm in col for gsm in gsm_ids)
        ]

        counts_df = counts_df[retained_cols]
        counts_df.index.name = None
        counts_df.insert(0, "Name", counts_df.index)
        counts_df.reset_index(drop=True, inplace=True)

        cache_path = geo_cache_dir / f"{accession}_{c.filename}_{id(c)}.parquet"
        modified_path = geo_cache_dir / f"{accession}_{c.filename}_{id(c)}_modified.parquet"
        counts_df.to_parquet(cache_path, compression="gzip")
        shutil.copy(cache_path, modified_path)

        meta_entry = {
            "path": cache_path,
            "modified_path": modified_path,
            "filename": c.filename,
            "normalization_type": c.normalization,
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

    df = df.drop("_probe_id", axis=1)
    print(f"Remapped across {len(results)} GPL platform(s)")
    return df



def fetch_microarray_matrix(meta: dict):
    dp_text = st.session_state["dp_text"]
    summ_norm = "unknown"
    if dp_text.strip():
        print("Classifying normalization type…")
        is_log, summ_norm = needs_log(dp_text)

    final_df = GLOBAL_GSE.pivot_samples(values="VALUE")[meta["gsm_ids"]]
    
    print(f"Microarray matrix contains {len(final_df.columns)} samples.")
    final_df = _annotate_matrix(final_df)
    final_df = reduce_matrix_memory(final_df)

    geo_cache_dir = Path("./geo_cache")
    geo_cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = geo_cache_dir / f"{gse_id}_microarray.parquet"
    modified_path = geo_cache_dir / f"{gse_id}_microarray_modified.parquet"

    final_df.to_parquet(cache_path, compression="gzip")
    shutil.copy(cache_path, modified_path)

    meta_entry = {
        "path": cache_path,
        "modified_path": modified_path,
        "filename": "microarray_matrix",
        "normalization_type": summ_norm,
        "n_genes": final_df.shape[0],
        "n_samples": final_df.shape[1] - 1,
    }
    del final_df
    gc.collect()
    return [meta_entry]

class StreamWrapper(io.RawIOBase):
    """Wraps a bytes-yielding generator into a read-only binary stream for Streamlit."""
    def __init__(self, generator):
        super().__init__()
        self.gen = generator
        self.leftover = b""

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        # Prevents crash when frameworks check seekability
        return False

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        # Ignore non-fatal seek attempts (e.g. seek(0, SEEK_CUR))
        if offset == 0 and whence in (io.SEEK_SET, io.SEEK_CUR):
            return 0
        raise io.UnsupportedOperation("stream is not seekable")

    def tell(self) -> int:
        return 0

    def readinto(self, b):
        try:
            while len(self.leftover) < len(b):
                chunk = next(self.gen)
                self.leftover += chunk
        except StopIteration:
            pass

        output, self.leftover = self.leftover[:len(b)], self.leftover[len(b):]
        b[:len(output)] = output
        return len(output)


def stream_parquet_to_tsv_gz(parquet_path: str):
    """Generator streaming TSV.GZ chunks directly from disk without stream corruption."""
    parquet_file = pq.ParquetFile(parquet_path)
    
    # Custom buffer to capture written gzip bytes incrementally
    class ChunkBuffer(io.RawIOBase):
        def __init__(self):
            self.chunks = []
        def writable(self):
            return True
        def write(self, b):
            self.chunks.append(bytes(b))
            return len(b)
        def get_and_clear(self):
            data = b"".join(self.chunks)
            self.chunks.clear()
            return data

    out_buf = ChunkBuffer()
    
    with gzip.GzipFile(fileobj=out_buf, mode="wb") as gz_file:
        for record_batch in parquet_file.iter_batches(batch_size=5000):
            df_chunk = record_batch.to_pandas()
            tsv_data = df_chunk.to_csv(sep="\t", index=False, header=False).encode('utf-8')
            gz_file.write(tsv_data)
            gz_file.flush()
            
            chunk_data = out_buf.get_and_clear()
            if chunk_data:
                yield chunk_data

    # Yield remaining bytes (including gzip trailer) after closing
    final_data = out_buf.get_and_clear()
    if final_data:
        yield final_data

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
            width='stretch',
            hide_index=True,
        )
        mortality = st.selectbox(
            "Event column",
            options=list(char_df.columns),
            key="survival_mortality_col"
        )

        # Iterate safely over unique values
        unique_vals = sorted(list(char_df[mortality].dropna().unique()))
        if [int(u) for u in unique_vals if str(u).isdigit()] == [0,1]:
            st.session_state.survival_df["death"] = char_df[mortality]
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

        st.session_state.survival_df[time_units] = char_df[t_mortality]

        st.write("Final time-to-event data: ")

        st.dataframe(
            st.session_state.survival_df.head(),
            width = 'stretch',
            hide_index = True
        )


        st.download_button(
            label="⬇ Download as .txt (tab-separated)",
            key = f"download_survival",
            data=st.session_state.survival_df.to_csv(sep="\t", index=False),
            file_name=f"{gse_id}_survival.txt",
            mime="text/plain",
        )


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

            modified_df = pq.read_table(modified_path)

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


            pq.write_table(
                modified_df, 
                modified_path, 
                compression="snappy"
            )

            st.session_state[f"selected_cols_dict_{gse_id}_{i}"] = {col: True for col in list(column_mapping.values())}
            st.session_state[f"selected_cols_dict_{gse_id}_{i}"]["Name"] = True
            st.rerun(scope="app")

@st.fragment
def column_selector(counts_path, i, gse_id):
    parquet_file = pq.ParquetFile(counts_path)
    all_columns = parquet_file.schema.names

    state_key = f"selected_cols_dict_{gse_id}_{i}"
    if state_key not in st.session_state:
        st.session_state[state_key] = {col: True for col in all_columns}

    search_term = st.text_input("Search columns", key=f"search_bar_{i}")

    filtered_columns = [
        col for col in all_columns
        if search_term.lower() in col.lower()
    ] if search_term else list(all_columns)

    def select_all():
        for col in filtered_columns:
            st.session_state[state_key][col] = True
            st.session_state[f"ui_{col}_{i}_key"] = True

    def deselect_all():
        for col in filtered_columns:
            st.session_state[state_key][col] = False
            st.session_state[f"ui_{col}_{i}_key"] = False

    col1, col2 = st.columns(2)
    with col1:
        st.button("Select All", on_click=select_all, key=f"select_all_{i}")
    with col2:
        st.button("Deselect All", on_click=deselect_all, key=f"deselect_all_{i}")

    with st.expander("Show All Columns", expanded=False):
        num_grid_cols = 4  
        grid_columns = st.columns(num_grid_cols)

        for idx, col in enumerate(filtered_columns):
            def update_single_col(c=col):
                st.session_state[state_key][c] = st.session_state[f"ui_{c}_{i}_key"]

            widget_key = f"ui_{col}_{i}_key"
            
            if col not in st.session_state[state_key]:
                st.session_state[state_key][col] = True
                
            if widget_key not in st.session_state:
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

def apply_norm(modified_path, norm, col_list, i, source_path):
    target_cols = set(col_list)

    schema = pq.read_schema(modified_path)
    base_schema = pq.read_schema(source_path)
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
        return pq.read_table(modified_path)

    numeric_df = pq.read_table(source_path, columns=base_gsm_cols)
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
        return pq.read_table(modified_path)

    del numeric_df
    gc.collect()

    df = pq.read_table(modified_path)

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

if st.button("Fetch & Build Matrix", type="primary"):
    st.session_state.run_pipeline = True
    st.session_state.result_lists = None

    for key in list(st.session_state.keys()):
        if str(key).startswith(("df_", "base_df_", "norm_type_", "preview_", "selected_cols_dict_", "pill_selector_")):
            del st.session_state[key]
    gc.collect()
    

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
        gsm_filter.extend(meta["gsm_gpl_dict"][g])

    survival_metadata_ui(char_df.loc[char_df.index.isin(gsm_filter) | ~char_df.index.str.startswith("GSM")])
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

        annotate_columns(i, char_df.loc[char_df.index.isin(gsm_filter) | ~char_df.index.str.startswith("GSM")], gse_id, meta_entry["path"], meta_entry["modified_path"])

        with st.expander("Change Normalization", expanded=False):
            st.write(st.session_state["dp_text"])
            is_log, summ_norm = needs_log(original_normalization)
            st.write("What normalization would you like to apply?")
            if original_normalization == "raw_counts":
                norm_suggestion = f"log2 or log2(cpm+1) because {summ_norm} data is raw."
            elif is_log:
                norm_suggestion = f"None because {summ_norm} already log scales the data."
            else:
                norm_suggestion = f"log2 because {summ_norm} is linearly scaled, potentially leading to a skewed distribution."

            st.radio(
                label=f"Suggestion: {norm_suggestion}",
                options=("none", "log2", "log10", "cpm", "log2(cpm+1)", "log10(cpm+1)"),
                key=f"norm_type_{i}"
            )

            st.write("Which columns would you like to apply it to?")
            st.session_state[f"selected_columns_{i}"] = column_selector(meta_entry["modified_path"], i, gse_id)

            submit_button = st.button(label="Renormalize", key=f"renormalize_{i}")

            if submit_button:
                modified_df = apply_norm(
                    meta_entry["modified_path"],
                    st.session_state[f"norm_type_{i}"],
                    st.session_state[f"selected_columns_{i}"],  
                    i,
                    meta_entry["path"]
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

        st.download_button(
            label="⬇ Download complete matrix (.txt.gz)",
            key=f"download_{i}",
            data=StreamWrapper(stream_parquet_to_tsv_gz(meta_entry["modified_path"])),
            file_name=f"{gse_id}_{original_normalization}{norm_suffix}.txt.gz",
            mime="application/gzip",
        )
