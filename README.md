# COMPASS Prep™ landing page

A desktop-first, responsive Streamlit landing page for **COMPASS Prep™ — data harmonization and preprocessing**.

The live public-data pathway links to [GEO-2-COMPASS](https://geo2compass.precsn.com/). RAW-2-COMPASS is presented as future downloadable standalone software, with its download control intentionally disabled and marked **Coming Soon**.

## Included

- `app.py` — the complete one-page Streamlit application
- `assets/compass-prep-background-4k.webp` — optimized 4K website background
- `assets/compass-prep-background-4k.png` — full-quality 4K background
- `assets/BACKGROUND_PROMPT.md` — the final background-generation prompt
- `requirements.txt` — Python dependency declaration
- `Procfile` — Railway-compatible process command
- `railway.json` — Railway build, start, health-check, and restart settings
- `runtime.txt` — Python runtime selection
- `.streamlit/config.toml` — Streamlit theme and server settings

## Run locally

Create and activate a Python virtual environment, then run:

```bash
pip install -r requirements.txt
streamlit run app.py
```

Streamlit will print a local address. Open it in a browser.

## Deploy on Railway

1. Create a new Git repository containing the files in this folder.
2. In Railway, choose **New Project → Deploy from GitHub repo** and select that repository.
3. Railway will detect the Python app and use `railway.json`/`Procfile` to start Streamlit on the assigned `$PORT`.
4. When deployment is healthy, open **Settings → Networking → Generate Domain**.

No secrets or environment variables are required for this landing page.

## Update later

When RAW-2-COMPASS is ready, replace the disabled button in `app.py` with a download link to the signed installer or release page. Keep installers outside the repository if they are large; a versioned release or object-storage URL is easier to maintain.

## Background asset

The background was generated specifically for this project using the built-in image-generation workflow, then upscaled to 3840×2160. The app embeds the optimized WebP locally, so it does not rely on an external image host.
