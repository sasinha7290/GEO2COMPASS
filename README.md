# GEO-2-COMPASS

GEO-2-COMPASS is a Streamlit application for retrieving gene-expression studies
from the NCBI Gene Expression Omnibus (GEO) and converting them into
analysis-ready expression matrices.

The application supports RNA-seq and microarray studies, sample metadata,
platform-based gene annotation, optional normalization, survival-metadata
generation, and compressed matrix downloads.

## Features

- Retrieve GEO Series records using a `GSE` accession
- Detect RNA-seq and microarray studies
- Select a platform when a study contains multiple GPL platforms
- Build expression matrices from GEO count files or supplementary files
- Map probe or gene identifiers to gene symbols when annotation is available
- Preview large matrices without sending the complete dataset to the browser
- Apply log, CPM, and log-CPM transformations
- Annotate sample columns using GEO characteristics
- Create survival-analysis metadata
- Download the complete matrix as a compressed, tab-separated file

## Requirements

- Python 3.11
- Internet access to NCBI GEO and g:Profiler
- At least 1 GB of memory is recommended for typical deployments

The application includes safety limits for hosted environments:

- Maximum individual download: 100 MB
- Maximum expression matrix: 20,000,000 cells
- Browser preview: 1,000 rows and 50 columns

These limits can be changed near the beginning of `app.py`, but increasing them
may require a larger Railway memory plan.

## Run locally

Clone or download the repository, then open a terminal in the project folder.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run app.py
```

On Windows PowerShell, activate the environment with:

```powershell
.\.venv\Scripts\Activate.ps1
```

Streamlit will display the local address in the terminal, normally
`http://localhost:8501`.

## Using GEO2COMPASS

1. Enter a GEO Series accession such as `GSE183620`.
2. Select **Load GEO metadata**.
3. If prompted, choose the GEO platform to process.
4. Select **Fetch & Build Matrix**.
5. Optionally annotate sample names or apply another normalization.
6. Download the complete matrix as a `.txt.gz` file.

Very large or unusually formatted GEO supplementary files may not be suitable
for automatic processing. The app reports a clear error when a file exceeds its
configured download or matrix limit.

## Deploy on Railway

1. Push every project file to a GitHub repository. Make sure the hidden
   `.streamlit` folder and the new `runtime.txt` file are included.
2. In Railway, create a project and select **Deploy from GitHub repo**.
3. Select this repository.
4. Railway will use `railway.toml` as the deployment configuration.
5. Wait for the `/_stcore/health` health check to pass before attaching a custom
   domain.

No user-defined environment variables are required for the standard
deployment. Railway supplies `PORT` automatically. Both `railway.toml` and the
`Procfile` bind Streamlit to:

```text
0.0.0.0:$PORT
```

The start command also sets:

```text
ARROW_DEFAULT_MEMORY_POOL=system
```

Do not replace `$PORT` with a fixed port in Railway.

## Railway configuration

The repository contains:

- `railway.toml` — start command, health check, and restart policy
- `Procfile` — compatible fallback start command
- `runtime.txt` and `.python-version` — Python 3.11 selection
- `.streamlit/config.toml` — headless server settings and disabled file watcher
- `requirements.txt` — pinned Python dependencies

If the deployment builds successfully but Railway reports that the application
failed to respond, inspect the deployment logs before changing DNS. A custom
domain cannot make an unhealthy application process respond.

## PyArrow crash protection

PyArrow 25.0.0 introduced a threaded memory-allocation crash that can be
triggered by Streamlit while serializing a pandas DataFrame. It appears as a
fatal Python segmentation fault in `pyarrow.pandas_compat`.

This repository applies both available protections:

- `pyarrow==24.0.0` is pinned in `requirements.txt`.
- `ARROW_DEFAULT_MEMORY_POOL=system` is set by the application and Railway
  start command.

Keep these protections until the affected PyArrow behavior is fixed and tested
with this application.

## Memory management

GEO studies can contain tens of millions of values. GEO2COMPASS reduces memory
use by:

- Processing one candidate matrix at a time
- Limiting concurrent downloads
- Streaming GEO metadata to a temporary file
- Reading TAR members individually
- Downcasting numeric columns when possible
- Caching only one parsed GEO accession
- Displaying a bounded preview instead of the complete matrix
- Producing gzip-compressed downloads

For datasets that exceed the built-in safety limit, use a Railway service with
more memory or process the dataset locally.

## Troubleshooting

### The app restarts with a segmentation fault

Confirm that the deployment installed `pyarrow==24.0.0` and that the Railway
start command still contains `ARROW_DEFAULT_MEMORY_POOL=system`. Redeploy
without the previous build cache if Railway reused an older dependency layer.

### Railway reports “Application failed to respond”

Check that Railway is using the start command in `railway.toml` and that the log
shows Streamlit listening on `0.0.0.0` with Railway's assigned port. Review the
earlier log entries for an application exception or out-of-memory termination.

### A matrix exceeds the safety limit

Run the application locally or increase the Railway service memory before
raising `MAX_MATRIX_CELLS` or `MAX_DOWNLOAD_BYTES` in `app.py`.

### Gene symbols are unavailable

Some GEO platforms do not provide a usable identifier-to-symbol mapping.
GEO2COMPASS retains the original identifiers when a reliable annotation cannot
be determined.

## Data sources

- [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/)
- [GEOparse](https://github.com/guma44/GEOparse)
- [g:Profiler](https://biit.cs.ut.ee/gprofiler/)

## License

No license file is currently included. Add a license before distributing or
accepting external contributions.

