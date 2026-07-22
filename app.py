from pathlib import Path
import tempfile
import GEOparse
import numpy as np
import pandas as pd
import streamlit as st
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup
import concurrent.futures
import gzip
from gprofiler import GProfiler 
import io
import json
import re
import tarfile
from dataclasses import dataclass, field
from typing import Optional
import numpy as np

st.set_page_config(page_title="GEO-2-COMPASS")
st.title("GEO-2-COMPASS")

gse_id = st.text_input("Enter GEO accession:", placeholder="e.g. GSE183620")

if not gse_id:
    st.stop()

gse_id = gse_id.strip().upper()

if "run_pipeline" not in st.session_state:
    st.session_state.run_pipeline = False
if "selected_columns" not in st.session_state:
    st.session_state.selected_columns = []
if "dp_text" not in st.session_state:
    st.session_state.dp_text = ""
if "char_df" not in st.session_state:
    st.session_state.char_df = pd.DataFrame()

@st.cache_data(show_spinner="Fetching GEO metadata…")
def fetch_metadata(gse_id: str) -> dict:
    try:
        #"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50081/soft/GSE50081_family.soft.gz"
        
        gpl_ids = list(GLOBAL_GSE.gpls.keys())
        gpl_titles = [
            GLOBAL_GSE.gpls[gpl_id].metadata.get('title', [''])[0] 
            for gpl_id in gpl_ids
        ]
        
        gsm_list = list(GLOBAL_GSE.gsms.keys())

        gsm_to_gpl = {}
        gsm_gpl_dict = {gpl_id: [] for gpl_id in gpl_ids}

        gsm_gpl_dict["Unknown"] = [] 
        
        for gsm_id, gsm_obj in GLOBAL_GSE.gsms.items():
            platform_list = gsm_obj.metadata.get('platform_id', [])
            associated_gpl = platform_list[0] if platform_list else ""
            gsm_to_gpl[gsm_id] = associated_gpl
            
            if associated_gpl in gsm_gpl_dict:
                gsm_gpl_dict[associated_gpl].append(gsm_id)
            elif associated_gpl:
                gsm_gpl_dict[associated_gpl] = [gsm_id]
            else:
                gsm_gpl_dict["Unknown"].append(gsm_id)
            
        gse_types = GLOBAL_GSE.metadata.get('type', [])
        gse_types_str = " ".join(gse_types).lower()
        study_type = "Microarray" if "array" in gse_types_str else "RNA-seq"
        
        taxon = ""
        if 'organism' in GLOBAL_GSE.metadata and GLOBAL_GSE.metadata['organism']:
            taxon = GLOBAL_GSE.metadata['organism'][0]
        elif GLOBAL_GSE.gsms:
            first_gsm = next(iter(GLOBAL_GSE.gsms.values()))
            organism_list = first_gsm.metadata.get('organism_ch1', [])
            if organism_list:
                taxon = organism_list[0]

        return {
            "title":      GLOBAL_GSE.metadata.get('title', [gse_id])[0],
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

def get_gse(gse_id):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/soft/{gse_id}_family.soft.gz"
    
    print(f" [Downloading] Fetching data directly from: {url}")
    
    try:
        # Stream the download
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        # Use a NamedTemporaryFile so it has a valid string path for GEOparse
        with tempfile.NamedTemporaryFile(suffix=".soft.gz", delete=True) as temp_file:
            # Write the downloaded chunks to the temporary file
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    temp_file.write(chunk)
            
            # Flush the buffer to ensure all data is written to disk
            temp_file.flush()
            
            print(f"Parsing data with GEOparse from temporary file...")
            # Pass the string path of the temporary file
            return GEOparse.get_GEO(filepath=temp_file.name, geotype = "GSE")
            
    except requests.exceptions.RequestException as e:
        st.error(f"Download failed: {e}")
        st.stop()
    except Exception as e:
        st.error(f"Failed to parse GEO data: {e}")
        st.stop()

def get_dp_and_char():
    first_gsm = GLOBAL_GSE.gsms[meta["gsm_ids"][0]]
    dp_text = first_gsm.metadata['data_processing'][0]
    st.session_state["dp_text"] = dp_text
    char_list = [x.split(": ")[0] for x in first_gsm.metadata["characteristics_ch1"]]
    char_df = pd.DataFrame(columns = char_list)
    for gsm in meta["gsm_ids"]:
        gsm_data = GLOBAL_GSE.gsms[gsm]
        new_row = {}
        for x in gsm_data.metadata["characteristics_ch1"]:
            if ": " in x:
                new_row[x.split(": ")[0]] = x.split(": ")[1]
            else:
                new_row[x.split(": ")[0]] = x.split(": ")[0]
        char_df.loc[gsm_data.metadata["geo_accession"][0]] = new_row
        if "treatment_protocol_ch1" in gsm_data.metadata.keys():
            char_df.loc[gsm_data.metadata["geo_accession"][0], "treatment_protocol"] = gsm_data.metadata["treatment_protocol_ch1"][0]
        if "growth_protocol_ch1" in gsm_data.metadata.keys():
            char_df.loc[gsm_data.metadata["geo_accession"][0], "growth_protocol"] = gsm_data.metadata["growth_protocol_ch1"][0]
    st.session_state["char_df"] = char_df

GLOBAL_GSE = get_gse(gse_id)
meta = fetch_metadata(gse_id)
get_dp_and_char()

if meta.get("error"):
    st.error(f"Failed to load metadata: {meta['error']}")
    st.stop()

st.subheader(meta["title"])
st.caption(
    f"Type: **{meta['type']}** · "
    f"Samples: **{len(meta['gsm_ids'])}** · "
    f"Taxon: {meta['taxon']}"
)

if len(meta["gpl_ids"]) > 1:
    selected_gpl = st.selectbox(
        "Multiple platforms detected. Which platform would you like to use?",
        tuple(f"{meta['gpl_ids'][i]}:  {meta['gpl_titles'][i]}" for i in range(len(meta['gpl_ids']))),
        index=None,
        placeholder="Select GPL"
    )
    if selected_gpl:
        selected_gpl = selected_gpl.split(":  ")[0]
        meta['gsm_ids'] = meta["gsm_gpl_dict"][selected_gpl]
        meta['gpl_ids'] = [selected_gpl]
else:
    selected_gpl = meta["gpl_ids"][0]

def needs_log(dp_text: str):
    yes_log = ["rma ", "(rma)", "lowess", "log", "vsn", "beadstudio","vsn","fhma","plier","quantile"]
    no_log = ["mas5", "mas 5"]
    for n in no_log:
        if n in dp_text.lower():
            return False, n
    for y in yes_log:
        if y in dp_text.lower(): 
            return True, y
    
    if "raw" in dp_text:
        return False, "raw"
    
    return False, dp_text

def to_cpm(df):
    numeric_df = df.apply(pd.to_numeric, errors="coerce")
    col_sums = numeric_df.sum(axis=0).replace(0, np.nan)
    return numeric_df.divide(col_sums, axis=1) * 1e6

# # # # # # # # # # # # # # # #
# GENE SYMBOL ANNOTATION  # # # 
# # # # # # # # # # # # # # # #

REJECT_PATTERNS = [
    re.compile(r'^(NM|NR|NP|XM|XR|XP|NG|NT|NC)_\d+(\.\d+)?$'),  
    re.compile(r'^AFFX-.*_at$'),
    re.compile(r'^.*_at$$'),
    re.compile(r'^[A-Z]{1,2}\d{6,}(\.\d+)?$'),
    re.compile(r'^GenMAPP'),                   
    re.compile(r'^ILMN_\d+$'),                                     
    re.compile(r'^ENSG\d+|^ENST\d+|^ENSMUSG\d+|^ENSMUST\d+|^ENSP\d+'),                  
    re.compile(r'^\d+(_[a-z]_at|_at|_x_at|_s_at)$'),             
    re.compile(r'^A_\d+_P\d+$'),                                   
    re.compile(r'^chr', re.IGNORECASE),                        
    re.compile(r'^\d+[pq]\d+'),                               
    re.compile(r'^GO:\d+$'),                                 
    re.compile(r'\.\d+$'),                                   
    re.compile(r'^\d'),
    re.compile(r'\s'),
    re.compile(r'\.'),                                   
    re.compile(r'_\d{4,}'),
    re.compile(r'scl\d+.'),
    re.compile(r'RefSeq'),
    re.compile(r'\+'),
    re.compile(r'^[agctAGCT]*$'),
    re.compile(r"^.{,3}$")                                       
]

ACCEPT_PATTERNS = [
    re.compile(r'^\d{7}[A-Za-z]\d{2}[Rr][Ii][Kk]$'),  
]

def is_gene_symbol(value):
    if "//" in str(value):
        value = str(value).split("//")
        v = value[1]
        v = str(v).strip()
        return (not any(p.search(v) for p in REJECT_PATTERNS)) or any(r.search(v) for r in ACCEPT_PATTERNS)

    else:
        v = str(value).strip()
        return (not any(p.search(v) for p in REJECT_PATTERNS)) or any(r.search(v) for r in ACCEPT_PATTERNS)
    
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

    if len(list(set(results))) < 100:
        return True

    col = list(set(col))
    return len(col) < 100

def gene_convert(gpl_df, gse, best_col, symbol_col):
    for gpl_name, gpl in getattr(gse, "gpls", {}).items():
        try:
            species = gpl.metadata.get('organism', [])[0]
        except (AttributeError, TypeError, KeyError):
            return gpl_df, symbol_col
    

    if (any(s in species.lower() for s in ["homo","sapiens"])):
        species = "hsapiens"
    elif (any(s in species.lower() for s in ["musculus","mus"])):
        species = "mmusculus"
    else:
        return gpl_df, symbol_col

    target_namespace = "HGNC" if species == "hsapiens" else "MGI"

    id = gpl_df[best_col].iloc[0]
    print(id.strip()[:3].lower() == "eg:")
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
        gpl_df["new_id_col"] = gpl_df[best_col].str.replace(r"_at$", "", regex=True)
        if (str(gpl_df[best_col].iloc[0]).strip()[:3].lower() == "eg:"):
            gpl_df["new_id_col"] = gpl_df[best_col].str.split(':').str[1]

        gp = GProfiler(return_dataframe=True)
        results = gp.convert(organism=species,
                    query=list(gpl_df["new_id_col"]),
                    numeric_namespace=numeric_namespace,
                    target_namespace=target_namespace)
        results_deduped = results.drop_duplicates(subset="incoming", keep="first")
        mapping = results_deduped.set_index("incoming")["name"]
        gpl_df["gene_symbol"] = gpl_df["new_id_col"].values 
        gpl_df["gene_symbol"] = [mapping.get(id_, float("nan")) for id_ in gpl_df["new_id_col"]]

        return gpl_df, "gene_symbol"
        
def get_gene_symbol_column(gse_id, counts_df, selected_gpl, probe_ids=None):
    gse = GLOBAL_GSE

    gpl_df = None
    for gpl_name, gpl in getattr(gse, "gpls", {}).items():
        if gpl_name == selected_gpl[0]:
            gpl_df = gpl.table
            print(gpl_name)
        #print(gpl.table)
        #print(gpl_df.empty)
    
    if probe_ids is not None:
        index_set = {_normalize_id(v) for v in probe_ids}
    else:
        index_set = {_normalize_id(v) for v in counts_df.index}

    if gpl_df.empty:
        gpl_df = pd.DataFrame({"id": list(index_set)})
        gpl_df["gene_symbol"] = gpl_df["id"]
        gpl_df, symbol_col = gene_convert(gpl_df, gse, "id", "gene_symbol")

        return gpl_df, "id", symbol_col
    
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
            if col_score > symbol_col_score:
                symbol_col_score = col_score
                symbol_col = col
        if symbol_col_score < 0.30: #arbitrary
            gpl_df, symbol_col = gene_convert(gpl_df, gse, best_col, symbol_col)
        
        return gpl_df, best_col, symbol_col
    else:
        raise ValueError(f"No platform data found for {gse}")
    

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
        df    = self.norm_df or self.counts_df
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

#bytes download

def _download_bytes(url: str) -> bytes:
    url = url.replace("ftp://", "https://", 1)
    head_resp = requests.head(url, timeout=30)
    head_resp.raise_for_status()
    
    file_size = int(head_resp.headers.get('Content-Length', 0))
    max_size = 100 * 1024 * 1024 
    
    if file_size > max_size:
        return None
    resp = requests.get(url, timeout=300)
    
    resp.raise_for_status()
    return resp.content

def _download_all(urls: list[str]) -> dict[str, bytes]:
    results: dict[str, bytes] = {}
    print(f"  Downloading {len(urls)} file(s) in parallel…")
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(64, len(urls))) as ex:
        fut_to_url = {ex.submit(_download_bytes, u): u for u in urls}
        for fut in concurrent.futures.as_completed(fut_to_url):
            url = fut_to_url[fut]
            fname = url.split("/")[-1]
            try:
                res = fut.result()
                results[url] = res
                if res is not None:
                    print(f"    ✓ {fname} ({len(results[url]) / 1e6:.1f} MB)")
                else:
                    print(f" {fname}: Exceeded size limit, returned None")
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


def _detect_sep(name: str) -> str:
    return "," if name.lower().endswith(".csv") else "\t"


def _unpack_tar(raw: bytes) -> dict[str, bytes]:
    buf = io.BytesIO(raw)
    members: dict[str, bytes] = {}
    with tarfile.open(fileobj=buf, mode="r:*") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            name = Path(member.name).name
            if _is_tabular(name):
                f = tf.extractfile(member)
                if f:
                    members[name] = f.read()
    return members

#tested
def _scrape_geo_download_page(accession):
    page_url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession}"
    print(f"  Scraping GEO download page: {page_url}")
 
    r = requests.get(page_url, timeout=60)
    r.raise_for_status()
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

def select_download_urls(candidates, accession, norm_type):
    submitter_candidates = []
    for c in candidates:
        if c.ncbi_data and "raw" in c.filename:
            return [c], "raw_counts"
        elif not c.ncbi_data:
            submitter_candidates.append(c)
    return [submitter_candidates[0]], norm_type

#constructing the dataframe

def _parse_tabular_series(raw: bytes, filename: str) -> tuple[pd.Series, str]:
    sep = _detect_sep(filename)
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

    s = pd.to_numeric(df[count_col], errors="coerce")
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
        members = _unpack_tar(raw)
        for mname, mbytes in members.items():
            
            mbytes, bare = _decompress(mbytes, mname)
            try:
                s, sname = _parse_tabular_series(mbytes, bare)
                s.name   = sname
                series_list.append(s)
            except Exception as e:
                print(f"    [warn] {mname}: {e}")
                return pd.DataFrame()
 
    if not series_list:
        return pd.DataFrame()
 
    print(f"\n  Merging {len(series_list)} sample series from TAR…")
    df = pd.concat(series_list, axis=1, join="outer")
    df = df.apply(pd.to_numeric, errors="coerce")
    df.index.name = "gene_id"
    print(df)
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
            sep = _detect_sep(filename)

            lower = filename.lower()
            if lower.endswith(".xls") or lower.endswith(".xlsx"):
                df = pd.read_excel(io.BytesIO(raw), index_col=0)
            else:
                df = pd.read_csv(
                    io.BytesIO(raw), sep=_detect_sep(filename),
                    index_col=0, engine="python", on_bad_lines="skip",
                )
            

            if df.columns.str.startswith('Unnamed').all():
                df.columns = df.iloc[0]

            df = df.apply(pd.to_numeric, errors = "coerce")
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

def _annotate_counts(accession, selected, counts_df, selected_gpl):
    if selected[0].ncbi_data:
        annot_url = "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
        url_bytes = _download_all([annot_url])
        raw = url_bytes[annot_url]
        filename = annot_url.split("file=")[-1]

        if filename.endswith(".gz"):
            raw = gzip.decompress(raw)
            filename = filename[:-3]

        annot_df = pd.read_csv(
            io.BytesIO(raw), sep=_detect_sep(filename),
            engine="python", on_bad_lines="skip",
        )
        
        id_col = "GeneID"
        symbol_col = "Symbol"
    else:
        try:
            annot_df, id_col, symbol_col = get_gene_symbol_column(accession, counts_df, selected_gpl)
        except:
            return counts_df

    mapping     = dict(zip(annot_df[id_col].astype(str), annot_df[symbol_col].astype(str)))

    mapping_lower = {str(k).lower(): v for k, v in mapping.items()}

    orig_index = counts_df.index.astype(str)
    new_index = orig_index.map(lambda x: mapping_lower.get(str(x).lower()))
    
    counts_df.index = new_index.where(
        new_index.notna() & (new_index != "nan"), orig_index
    )
    counts_df.index.name = "gene_id"
    print(f"  ✓ Index remapped using '{symbol_col}'")
    return counts_df


def _read_dp_text(geo, accession):
    if accession.startswith("GSE"):
        gsm_dict = geo.gsms
    else:
        gsm_dict = {accession: geo}
 
    first_gsm = next(iter(gsm_dict.values()))
    return "\n".join(first_gsm.metadata.get("data_processing", []))

def fetch_and_normalize(
    selected_gpl,
    accession:     str,
    geo_cache_dir: str | Path = "./geo_cache",
    save_output:   bool       = False,
    output_dir:    str | Path = "./geo_output",
) -> list[GseResult]:

    geo_cache_dir = Path(geo_cache_dir)
    geo_cache_dir.mkdir(parents=True, exist_ok=True)
    accession = accession.strip().upper()
 
    print(f"\n{'═'*60}")
    print(f"  GEO accession: {accession}")
    print(f"{'═'*60}")
    dp_text  = st.session_state["dp_text"]
 
    candidates = _scrape_geo_download_page(accession)
    if not candidates:
        return [GseResult(
            accession=accession, normalization_type="unknown",
            effective_norm="unknown"
        )]
    
    num_supplementary = 0

    for c in candidates:
        if not c.ncbi_data:
            num_supplementary += 1

    results = []
    
    for c in candidates:
        selected = [c]
        print(f"\n  Selected {len(selected)} file(s) to download:")
        for m in selected:
            print(f"    {'[TAR]' if m.is_tar else '     '} {m.filename}")
    
        url_bytes = _download_all([m.url for m in selected])

        if url_bytes.get(c.url) is None:
            print(f"Exceeded file size, skipping this file {c.url}")
            continue

        if num_supplementary == 1 and not c.ncbi_data:
            c.normalization = classify_normalization(dp_text)

        print("\n  Building count matrix…")
        counts_df = _fetch_counts_df(accession, selected, url_bytes)
        
        if counts_df.empty:
            continue
        
        counts_df = _annotate_counts(accession, selected, counts_df, selected_gpl)
        
        result = GseResult(
            accession=accession,
            normalization_type=c.normalization,
            effective_norm=None,
            counts_df=counts_df,
            norm_df=None
        )

        results.append(result)

    if save_output:
        results.save(output_dir)

    return results

def fetch_rnaseq_matrix(gse_id: str, meta):
    try:
        results = fetch_and_normalize(
            selected_gpl=meta["gpl_ids"],
            accession=gse_id,
            geo_cache_dir="./geo_cache",
            save_output=False,
        )
    except Exception as e:
        st.error(f"geo_rnaseq_normalizer error: {e}")
        return None
    
    result_lists = []

    for result in results:
        filtered_cols = [col for col in list(result.counts_df.columns) if any(gsm in col for gsm in meta["gsm_ids"]) and "gsm" in col.lower() or "gsm" not in col.lower()]
        result.counts_df = result.counts_df[filtered_cols]
        result.norm_df = result.counts_df

        result.effective_norm = "none"

        df = result.norm_df.copy()
        df.index.name = None
        df.insert(0, "Name", df.index)
        df = df.reset_index(drop=True)

        st.info(
            f"Original normalization: **{result.normalization_type}** "
        )

        #print(df)

        result_list = [df, result.normalization_type, result.effective_norm]
        result_lists.append(result_list)

    return result_lists

# # # # # # # # # # # # # # # #
# MICROARRAY DATASET HANDLING # 
# # # # # # # # # # # # # # # #

def _annotate_matrix(df: pd.DataFrame) -> pd.DataFrame:
    df = df.reset_index()
    df.rename(columns={df.columns[0]: "_probe_id"}, inplace=True)

    probe_ids = [str(pid).strip() for pid in df["_probe_id"]]
    looks_like_symbols = sum([int(is_gene_symbol(pid)) for pid in probe_ids[:min(100, len(probe_ids))]])

    if looks_like_symbols:
        df.insert(0, "Name", df["_probe_id"])
        df = df.drop("_probe_id", axis=1)
        return df

    NULL_VALUES = {"---", "na", "", "null", "nan"}

    gpl_table, matching_col, symbol_col = get_gene_symbol_column(
        gse_id, df, meta["gpl_ids"], probe_ids=df["_probe_id"]
    )
    print(f"  → matching_col='{matching_col}'  symbol_col='{symbol_col}'")

    master_mapping: dict[str, str] = {}
    if symbol_col in gpl_table.columns and matching_col in gpl_table.columns:
        for probe, sym in zip(
            gpl_table[matching_col].astype(str).str.strip(),
            gpl_table[symbol_col].astype(str).str.strip(),
        ):
            if probe and sym and sym.lower() not in NULL_VALUES:
                if "//" in sym:
                    sym = re.split(r'///|//', sym)[1].strip()
                master_mapping[probe] = sym

    df["Name"] = df["_probe_id"].astype(str).str.strip().map(master_mapping)

    before = len(df)
    df = df.dropna(subset=["Name"])
    after = len(df)
    if before != after:
        print(f"dropped {before - after} unmapped rows ({after} remaining)")

    df = df.drop("_probe_id", axis=1)
    print(f"Remapped using '{matching_col}' → '{symbol_col}'")
    return df

def fetch_microarray_matrix(meta: dict):
    dp_text = st.session_state["dp_text"]
    if dp_text.strip():
        print("Classifying normalization type…")
        is_log, summ_norm = needs_log(dp_text)

    final_df = GLOBAL_GSE.pivot_samples(values="VALUE")[meta["gsm_ids"]]
    print("\t \t COLUMNS \t \t")
    print(final_df.columns)
    print(final_df)
    print(f"\n Number of Columns: \t \t {len(final_df.columns)}")
    final_df = _annotate_matrix(final_df)

    norm_type = "none"

    return [[final_df, summ_norm, norm_type]]

# # # # # # # # # # # # # # # #
# USER INTERFACE  # # # # # # # 
# # # # # # # # # # # # # # # #

@st.fragment
def survival_metadata_ui(char_df):
    if "survival_df" not in st.session_state:
        st.session_state.survival_df = pd.DataFrame(index=char_df.index)
        st.session_state.survival_df["GSM"] = char_df.index

    with st.expander("Generate Survival Metadata File"):
        st.dataframe(
            char_df.head(),
            width='stretch',
            hide_index=True,
        )
        mortality = st.selectbox(
            "Mortality column",
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
                    label=f"Label '{i}' as alive (0) or dead (1)",
                    options=(0, 1),
                    key=f"alive_dead_{i}"
                )
            st.session_state.survival_df["death"] = char_df[mortality].map(mapping)

        t_mortality = st.selectbox(
            "Time to mortality column",
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

        st.write("Final survival metadata file: ")

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

def column_selector(counts_df, i, gse_id):
    state_key = f"selected_cols_dict_{gse_id}_{i}"
    if state_key not in st.session_state:
        st.session_state[state_key] = {col: True for col in counts_df.columns}

    search_term = st.text_input("Search columns", key=f"search_bar_{i}")

    filtered_columns = [
        col for col in counts_df.columns
        if search_term.lower() in col.lower()
    ] if search_term else list(counts_df.columns)

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
        col for col in counts_df.columns
        if st.session_state[state_key].get(col, True)
    ]
    
    return selected_columns

def apply_norm(df, norm, col_list, i):
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    numeric_cols = [n for n in numeric_cols if n in col_list]
    
    base_gsm_cols = []
    valid_numeric_cols = []
    
    for col in numeric_cols:
        gsm_id = col.split("_")[-1]  
        
        if gsm_id in st.session_state[f"base_df_{i}"].columns:
            base_gsm_cols.append(gsm_id)
            valid_numeric_cols.append(col)
        elif col in st.session_state[f"base_df_{i}"].columns:
            base_gsm_cols.append(col)
            valid_numeric_cols.append(col)

    if not base_gsm_cols:
        return df  
    numeric_df = st.session_state[f"base_df_{i}"][base_gsm_cols].copy()
    
    numeric_df.columns = valid_numeric_cols
    if norm == "none":
        return df.assign(**{
            col: numeric_df[col] for col in valid_numeric_cols
        })

    if norm == "log2":
        return df.assign(**{
            col: np.log2(numeric_df[col] + 1) for col in valid_numeric_cols
        })

    if norm == "log10":
        return df.assign(**{
            col: np.log10(numeric_df[col] + 1) for col in valid_numeric_cols
        })

    if norm == "cpm":
        cpm = numeric_df.div(numeric_df.sum(axis=0), axis=1) * 1e6
        return df.assign(**cpm.to_dict("series"))

    if norm == "log2(cpm+1)":
        cpm = numeric_df.div(numeric_df.sum(axis=0), axis=1) * 1e6
        return df.assign(**{
            col: np.log2(cpm[col] + 1) for col in valid_numeric_cols
        })

    if norm == "log10(cpm+1)":
        cpm = numeric_df.div(numeric_df.sum(axis=0), axis=1) * 1e6
        return df.assign(**{
            col: np.log10(cpm[col] + 1) for col in valid_numeric_cols
        })
        
    return df
    
if "result_lists" not in st.session_state:
    st.session_state.result_lists = None

if st.button("Fetch & Build Matrix", type="primary"):
    st.session_state.run_pipeline = True
    st.session_state.result_lists = None

    for key in list(st.session_state.keys()):
        if key.startswith("df_") or key.startswith("base_df_") or key.startswith("norm_type_"):
            del st.session_state[key]
    

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

        if result_lists is None:
            status.update(label="Pipeline failed.", state="error")
            st.stop()
        st.session_state.result_lists = result_lists
        status.update(label="Done!", state="complete")

if st.session_state.result_lists is not None:
    result_lists = st.session_state.result_lists
    survival_metadata_ui(st.session_state.char_df) 

    for i in range(len(result_lists)):
        if f"norm_type_{i}" not in st.session_state:
            st.session_state[f"norm_type_{i}"] = "none"
        result_list = result_lists[i]
        combined_df = result_list[0]
        if f"base_df_{i}" not in st.session_state:
            st.session_state[f"base_df_{i}"] = combined_df
    
        df_key = f"df_{i}"

        if df_key not in st.session_state:
            combined_df_copy = combined_df.copy()
            combined_df_copy.insert(0, 'Name', combined_df_copy.pop('Name'))
            st.session_state[df_key] = combined_df_copy

        df = st.session_state[df_key]
        
        original_normalization = result_list[1]
        effective_normalization = result_list[2]
        n_genes   = combined_df.shape[0]
        n_samples = combined_df.shape[1] - (1 if "Name" in combined_df.columns else 0)
        st.markdown(body = "<hr>", unsafe_allow_html= True)
        c4, c5 = st.columns(2)
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
        
        with c4:
            st.markdown(body=f'<div class = "container"><p class = "normalization">Original Normalization</p> <p class = "normtext">{original_normalization}</p></div>', unsafe_allow_html = True)

        with c5:
            st.markdown(body=f'<div class = "container"><p class = "normalization">Applied Normalization</p> <p class = "normtext">{st.session_state[f"norm_type_{i}"].upper()}</p></div>', unsafe_allow_html = True)

        c1, c2, c3 = st.columns(3)
        c1.metric("Genes",   f"{n_genes:,}")
        c2.metric("Samples", n_samples)

        with st.expander("Annotate Columns", expanded = False):
            char_df = st.session_state.char_df
            #print(char_df)
            st.write("Replace GSM column names with sample data")
            
            st.pills("Annotate by:", char_df.columns, selection_mode="multi", key = f"pill_selector_{i}")
            #print(st.session_state[f"pill_selector_{i}"])
            current_selection = st.session_state.get(f"pill_selector_{i}") or []
            if current_selection:
                example_name = "_".join([str(char_df[x].iloc[0]) for x in st.session_state[f"pill_selector_{i}"]])
                st.write(f"Example sample name: {example_name}_{char_df.index[0]}")
            else:
                st.write(char_df.index[0])
            if st.button("Annotate Columns", key = f"annotate_columns_{i}"):
                column_mapping = {}
                if not current_selection:
                    for gsm in char_df.index:
                        column_mapping[gsm] = str(gsm)
                else:
                    for gsm in char_df.index:
                        column_mapping[gsm] = "_".join([str(char_df[x].loc[gsm]) for x in current_selection]) + f"_{gsm}"
                
                base_df = st.session_state[f"base_df_{i}"].copy()
                
                current_norm = st.session_state.get(f"norm_type_{i}", "none")
                if current_norm != "none":
                    base_df = apply_norm(base_df, current_norm, st.session_state.get("selected_columns", list(base_df.columns)), i)
                    
                
                st.session_state[df_key] = base_df.rename(columns=column_mapping)

                st.session_state[f"selected_cols_dict_{gse_id}_{i}"] = {col: True for col in list(column_mapping.values())}
                st.session_state[f"selected_cols_dict_{gse_id}_{i}"]["Name"] = True
                st.rerun()
        
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
                label = f"Suggestion: {norm_suggestion}",
                options = ("none", "log2", "log10", "cpm", "log2(cpm+1)", "log10(cpm+1)"),
                key = f"norm_type_{i}"
            )

            st.write("Which columns would you like to apply it to?")
            st.session_state["selected_columns"] = column_selector(df,i,gse_id)
                
            submit_button = st.button(label="Renormalize", key = f"renormalize_{i}")

            if submit_button:
                st.session_state[df_key] = apply_norm(
                    st.session_state[df_key], 
                    st.session_state[f"norm_type_{i}"], 
                    st.session_state.selected_columns, 
                    i
                )
                df = st.session_state[df_key]
        
        st.dataframe(
            df,
            width='stretch',
            hide_index=True,
            column_config={
                "Name": st.column_config.TextColumn("Name", pinned=True, width="small"),
            },
        )

        tsv = df.to_csv(sep="\t", index=False)
        norm_suffix = "_" + st.session_state[f"norm_type_{i}"] if st.session_state[f"norm_type_{i}"] != "none" else ""
        st.download_button(
            label="⬇ Download as .txt (tab-separated)",
            key = f"download_{i}",
            data=tsv,
            file_name=f"{gse_id}_{original_normalization}{norm_suffix}.txt",
            mime="text/plain",
        )
