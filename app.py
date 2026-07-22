"""COMPASS Prep landing page.

Run locally with:
    streamlit run app.py
"""

from __future__ import annotations

import base64
from pathlib import Path

import streamlit as st


APP_DIR = Path(__file__).resolve().parent
BACKGROUND_PATH = APP_DIR / "assets" / "compass-prep-background-4k.webp"


def asset_data_uri(path: Path, mime_type: str) -> str:
    """Return a local image as an embeddable data URI."""
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime_type};base64,{encoded}"


st.set_page_config(
    page_title="COMPASS Prep | Coding-free preprocessing",
    page_icon="🧭",
    layout="wide",
    initial_sidebar_state="collapsed",
)

background_uri = asset_data_uri(BACKGROUND_PATH, "image/webp")

st.html(
    f"""
    <style>
        :root {{
            --navy: #08264a;
            --blue: #1677c8;
            --teal: #159f9a;
            --gold: #d89a2b;
            --ink: #10243d;
            --muted: #5f7185;
            --line: rgba(10, 73, 123, 0.14);
            --ice: #f4fbfc;
            --white: #ffffff;
        }}

        * {{ box-sizing: border-box; }}

        html {{ scroll-behavior: smooth; }}

        body {{
            margin: 0;
            background: #f7fbfc;
            color: var(--ink);
        }}

        body, button, a {{
            font-family: Inter, ui-sans-serif, -apple-system, BlinkMacSystemFont,
                         "Segoe UI", sans-serif;
        }}

        [data-testid="stHeader"],
        [data-testid="stToolbar"],
        [data-testid="stDecoration"],
        [data-testid="stStatusWidget"],
        footer {{
            display: none !important;
        }}

        [data-testid="stAppViewContainer"] {{
            background: #f7fbfc;
        }}

        [data-testid="stMain"] {{
            padding: 0;
        }}

        .block-container {{
            max-width: none;
            padding: 0 !important;
        }}

        .compass-page {{
            min-height: 100vh;
            overflow: hidden;
            background: linear-gradient(180deg, #ffffff 0%, #f5fbfc 62%, #ffffff 100%);
        }}

        .site-nav {{
            position: absolute;
            z-index: 10;
            top: 0;
            left: 50%;
            width: min(1180px, calc(100% - 64px));
            height: 84px;
            transform: translateX(-50%);
            display: flex;
            align-items: center;
            justify-content: space-between;
            border-bottom: 1px solid rgba(8, 38, 74, 0.09);
        }}

        .brand {{
            display: flex;
            align-items: center;
            gap: 12px;
            color: var(--navy);
            font-size: 17px;
            font-weight: 760;
            letter-spacing: -0.02em;
            text-decoration: none;
        }}

        .brand-mark {{
            position: relative;
            width: 36px;
            height: 36px;
            border: 2px solid rgba(22, 119, 200, 0.45);
            border-radius: 50%;
            background: rgba(255, 255, 255, 0.82);
            box-shadow: 0 5px 18px rgba(22, 119, 200, 0.11);
        }}

        .brand-mark::before {{
            content: "";
            position: absolute;
            top: 5px;
            left: 14px;
            width: 7px;
            height: 22px;
            border-radius: 8px 8px 2px 2px;
            background: linear-gradient(180deg, var(--gold) 0 48%, var(--blue) 48% 100%);
            transform: rotate(38deg);
            transform-origin: center;
        }}

        .brand-mark::after {{
            content: "";
            position: absolute;
            inset: 13px;
            border-radius: 50%;
            background: white;
        }}

        .nav-links {{
            display: flex;
            align-items: center;
            gap: 30px;
        }}

        .nav-links a {{
            color: #38516a;
            font-size: 14px;
            font-weight: 630;
            text-decoration: none;
            transition: color .2s ease;
        }}

        .nav-links a:hover {{ color: var(--blue); }}

        .hero {{
            position: relative;
            min-height: 680px;
            padding: 146px max(32px, calc((100vw - 1180px) / 2)) 76px;
            background-image:
                linear-gradient(90deg,
                    rgba(255,255,255,.98) 0%,
                    rgba(255,255,255,.94) 37%,
                    rgba(255,255,255,.53) 66%,
                    rgba(255,255,255,.08) 100%),
                url("{background_uri}");
            background-size: cover;
            background-position: center center;
            border-bottom: 1px solid var(--line);
        }}

        .hero::after {{
            content: "";
            position: absolute;
            left: 0;
            right: 0;
            bottom: 0;
            height: 120px;
            background: linear-gradient(180deg, rgba(247,251,252,0), #f7fbfc);
            pointer-events: none;
        }}

        .hero-copy {{
            position: relative;
            z-index: 2;
            width: min(650px, 58vw);
        }}

        .eyebrow {{
            display: inline-flex;
            align-items: center;
            gap: 9px;
            margin: 0 0 22px;
            color: var(--blue);
            font-size: 12px;
            font-weight: 800;
            letter-spacing: .16em;
            text-transform: uppercase;
        }}

        .eyebrow::before {{
            content: "";
            width: 24px;
            height: 2px;
            background: linear-gradient(90deg, var(--teal), var(--gold));
        }}

        .hero h1 {{
            margin: 0;
            color: var(--navy);
            font-size: clamp(52px, 5.6vw, 84px);
            line-height: .98;
            letter-spacing: -0.058em;
        }}

        .hero h1 sup {{
            position: relative;
            top: -1.8em;
            margin-left: 4px;
            color: var(--blue);
            font-size: 15px;
            letter-spacing: 0;
        }}

        .hero-tagline {{
            margin: 25px 0 0;
            color: #24445f;
            font-size: clamp(20px, 2vw, 27px);
            font-weight: 540;
            line-height: 1.35;
            letter-spacing: -0.025em;
        }}

        .hero-body {{
            max-width: 600px;
            margin: 18px 0 0;
            color: var(--muted);
            font-size: 16px;
            line-height: 1.72;
        }}

        .hero-actions {{
            display: flex;
            align-items: center;
            gap: 22px;
            margin-top: 32px;
        }}

        .primary-link {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            min-height: 50px;
            padding: 0 24px;
            border-radius: 9px;
            background: linear-gradient(135deg, var(--blue), #0d5fa6);
            box-shadow: 0 12px 28px rgba(22, 119, 200, .23);
            color: white !important;
            font-size: 14px;
            font-weight: 740;
            text-decoration: none;
            transition: transform .2s ease, box-shadow .2s ease;
        }}

        .primary-link:hover {{
            transform: translateY(-2px);
            box-shadow: 0 16px 34px rgba(22, 119, 200, .3);
        }}

        .primary-link span {{
            margin-left: 10px;
            font-size: 19px;
            line-height: 1;
        }}

        .quiet-link {{
            color: #3f607d !important;
            font-size: 14px;
            font-weight: 670;
            text-decoration: none;
            border-bottom: 1px solid rgba(63, 96, 125, .34);
        }}

        .trust-row {{
            display: flex;
            flex-wrap: wrap;
            gap: 12px 22px;
            margin-top: 36px;
            color: #516b83;
            font-size: 12px;
            font-weight: 650;
        }}

        .trust-row span::before {{
            content: "✓";
            margin-right: 7px;
            color: var(--teal);
            font-weight: 900;
        }}

        .section {{
            width: min(1180px, calc(100% - 64px));
            margin: 0 auto;
            padding: 84px 0;
        }}

        .section-kicker {{
            margin: 0 0 10px;
            color: var(--teal);
            font-size: 12px;
            font-weight: 800;
            letter-spacing: .15em;
            text-transform: uppercase;
        }}

        .section-title {{
            max-width: 720px;
            margin: 0;
            color: var(--navy);
            font-size: clamp(31px, 3vw, 44px);
            line-height: 1.12;
            letter-spacing: -.04em;
        }}

        .section-lede {{
            max-width: 660px;
            margin: 16px 0 0;
            color: var(--muted);
            font-size: 16px;
            line-height: 1.65;
        }}

        .path-grid {{
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 24px;
            margin-top: 40px;
        }}

        .path-card {{
            position: relative;
            min-height: 355px;
            padding: 34px;
            overflow: hidden;
            border: 1px solid var(--line);
            border-radius: 20px;
            background: rgba(255,255,255,.92);
            box-shadow: 0 20px 60px rgba(13, 66, 105, .08);
        }}

        .path-card.live {{
            background:
                radial-gradient(circle at 93% 8%, rgba(32,167,160,.16), transparent 32%),
                white;
        }}

        .path-card.soon {{
            background:
                radial-gradient(circle at 94% 8%, rgba(216,154,43,.14), transparent 31%),
                linear-gradient(135deg, #ffffff, #fbfcfc);
        }}

        .card-status {{
            display: inline-flex;
            align-items: center;
            gap: 8px;
            margin-bottom: 28px;
            padding: 7px 11px;
            border-radius: 100px;
            color: #107d78;
            background: #e8f8f6;
            font-size: 11px;
            font-weight: 820;
            letter-spacing: .08em;
            text-transform: uppercase;
        }}

        .card-status::before {{
            content: "";
            width: 7px;
            height: 7px;
            border-radius: 50%;
            background: var(--teal);
            box-shadow: 0 0 0 4px rgba(21,159,154,.12);
        }}

        .card-status.coming {{
            color: #8b621a;
            background: #fcf4e5;
        }}

        .card-status.coming::before {{
            background: var(--gold);
            box-shadow: 0 0 0 4px rgba(216,154,43,.12);
        }}

        .path-card h3 {{
            margin: 0;
            color: var(--navy);
            font-size: 28px;
            line-height: 1.2;
            letter-spacing: -.035em;
        }}

        .path-label {{
            margin: 7px 0 0;
            color: var(--blue);
            font-size: 12px;
            font-weight: 800;
            letter-spacing: .11em;
            text-transform: uppercase;
        }}

        .path-card p {{
            margin: 18px 0 22px;
            color: var(--muted);
            font-size: 15px;
            line-height: 1.62;
        }}

        .mini-list {{
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-bottom: 26px;
        }}

        .mini-list span {{
            padding: 7px 10px;
            border: 1px solid rgba(13, 95, 166, .12);
            border-radius: 7px;
            color: #49637b;
            background: #f8fbfd;
            font-size: 11px;
            font-weight: 650;
        }}

        .card-action {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            min-height: 46px;
            padding: 0 18px;
            border: 0;
            border-radius: 8px;
            color: white !important;
            background: var(--navy);
            font-size: 13px;
            font-weight: 750;
            text-decoration: none;
            transition: background .2s ease, transform .2s ease;
        }}

        .card-action:hover {{
            transform: translateY(-1px);
            background: var(--blue);
        }}

        .disabled-action {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            min-height: 46px;
            padding: 0 18px;
            border: 1px solid #dbe2e7;
            border-radius: 8px;
            color: #8795a2;
            background: #edf1f3;
            font-size: 13px;
            font-weight: 750;
            cursor: not-allowed;
        }}

        .disabled-action small {{
            margin-left: 9px;
            padding-left: 9px;
            border-left: 1px solid #cfd8de;
            font-size: 10px;
            letter-spacing: .07em;
            text-transform: uppercase;
        }}

        .workflow-wrap {{
            padding: 78px max(32px, calc((100vw - 1180px) / 2));
            background: var(--navy);
            color: white;
        }}

        .workflow-top {{
            display: flex;
            align-items: end;
            justify-content: space-between;
            gap: 40px;
        }}

        .workflow-top h2 {{
            max-width: 600px;
            margin: 0;
            font-size: clamp(30px, 3vw, 43px);
            line-height: 1.13;
            letter-spacing: -.04em;
        }}

        .workflow-top p {{
            max-width: 420px;
            margin: 0;
            color: #b9cedd;
            font-size: 15px;
            line-height: 1.6;
        }}

        .workflow {{
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 0;
            margin-top: 46px;
            border: 1px solid rgba(255,255,255,.14);
            border-radius: 16px;
            overflow: hidden;
        }}

        .step {{
            position: relative;
            min-height: 180px;
            padding: 28px;
            border-right: 1px solid rgba(255,255,255,.14);
            background: rgba(255,255,255,.035);
        }}

        .step:last-child {{ border-right: 0; }}

        .step-number {{
            color: #5bcac4;
            font-size: 12px;
            font-weight: 850;
            letter-spacing: .1em;
        }}

        .step h3 {{
            margin: 18px 0 8px;
            color: white;
            font-size: 19px;
        }}

        .step p {{
            margin: 0;
            color: #aac1d2;
            font-size: 14px;
            line-height: 1.55;
        }}

        .principles {{
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 24px;
            margin-top: 44px;
        }}

        .principle {{
            padding-top: 20px;
            border-top: 2px solid #dcebef;
        }}

        .principle strong {{
            display: block;
            color: var(--navy);
            font-size: 16px;
        }}

        .principle span {{
            display: block;
            margin-top: 8px;
            color: var(--muted);
            font-size: 13px;
            line-height: 1.5;
        }}

        .site-footer {{
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 24px;
            width: min(1180px, calc(100% - 64px));
            margin: 0 auto;
            padding: 28px 0 34px;
            border-top: 1px solid var(--line);
            color: #718294;
            font-size: 12px;
        }}

        .footer-brand {{
            color: var(--navy);
            font-weight: 760;
        }}

        .footer-links {{
            display: flex;
            gap: 22px;
        }}

        .footer-links a {{
            color: #566f87;
            text-decoration: none;
        }}

        .footer-links a:hover {{ color: var(--blue); }}

        @media (max-width: 900px) {{
            .site-nav {{ width: calc(100% - 40px); }}
            .nav-links a:not(:last-child) {{ display: none; }}
            .hero {{
                min-height: 660px;
                padding: 132px 24px 62px;
                background-position: 62% center;
            }}
            .hero-copy {{ width: min(610px, 88vw); }}
            .path-grid {{ grid-template-columns: 1fr; }}
            .workflow-top {{ align-items: start; flex-direction: column; }}
            .principles {{ grid-template-columns: repeat(2, 1fr); }}
        }}

        @media (max-width: 620px) {{
            .site-nav {{ height: 70px; width: calc(100% - 32px); }}
            .brand {{ font-size: 15px; }}
            .brand-mark {{ width: 32px; height: 32px; }}
            .brand-mark::before {{ top: 4px; left: 12px; }}
            .brand-mark::after {{ inset: 11px; }}
            .hero {{
                min-height: 700px;
                padding: 116px 20px 58px;
                background-image:
                    linear-gradient(180deg, rgba(255,255,255,.98) 0%, rgba(255,255,255,.92) 55%, rgba(255,255,255,.75) 100%),
                    url("{background_uri}");
            }}
            .hero-copy {{ width: 100%; }}
            .hero h1 {{ font-size: 50px; }}
            .hero h1 sup {{ top: -1.5em; }}
            .hero-actions {{ align-items: flex-start; flex-direction: column; gap: 16px; }}
            .section {{ width: calc(100% - 40px); padding: 62px 0; }}
            .path-card {{ padding: 27px; }}
            .workflow-wrap {{ padding: 62px 20px; }}
            .workflow {{ grid-template-columns: 1fr; }}
            .step {{ min-height: auto; border-right: 0; border-bottom: 1px solid rgba(255,255,255,.14); }}
            .step:last-child {{ border-bottom: 0; }}
            .principles {{ grid-template-columns: 1fr; }}
            .site-footer {{ align-items: flex-start; flex-direction: column; width: calc(100% - 40px); }}
        }}

        @media (prefers-reduced-motion: reduce) {{
            html {{ scroll-behavior: auto; }}
            * {{ transition: none !important; }}
        }}
    </style>

    <main class="compass-page" id="top">
        <nav class="site-nav" aria-label="Primary navigation">
            <a class="brand" href="#top" aria-label="COMPASS Prep home">
                <span class="brand-mark" aria-hidden="true"></span>
                <span>COMPASS Prep</span>
            </a>
            <div class="nav-links">
                <a href="#pathways">Pathways</a>
                <a href="#workflow">How it works</a>
                <a href="https://geo2compass.precsn.com/" target="_blank" rel="noopener noreferrer">Launch GEO-2-COMPASS ↗</a>
            </div>
        </nav>

        <header class="hero">
            <div class="hero-copy">
                <p class="eyebrow">Data harmonization &amp; preprocessing</p>
                <h1>COMPASS Prep<sup>™</sup></h1>
                <p class="hero-tagline">Coding-free preprocessing for deterministic biology.</p>
                <p class="hero-body">
                    Transform public transcriptomic datasets—or, soon, your own raw sequencing outputs—into standardized, COMPASS-ready inputs through clear and reproducible workflows.
                </p>
                <div class="hero-actions">
                    <a class="primary-link" href="https://geo2compass.precsn.com/" target="_blank" rel="noopener noreferrer">
                        Start with public data <span aria-hidden="true">→</span>
                    </a>
                    <a class="quiet-link" href="#pathways">Explore both pathways</a>
                </div>
                <div class="trust-row" aria-label="Platform benefits">
                    <span>Coding-free</span>
                    <span>Standardized outputs</span>
                    <span>Reproducible workflow</span>
                </div>
            </div>
        </header>

        <section class="section" id="pathways">
            <p class="section-kicker">Choose your starting point</p>
            <h2 class="section-title">One preparation layer. Two routes into COMPASS.</h2>
            <p class="section-lede">
                Use the live web workflow for public GEO studies today. A downloadable application for processing custom datasets locally is in development.
            </p>

            <div class="path-grid">
                <article class="path-card live">
                    <div class="card-status">Available now</div>
                    <h3>GEO-2-COMPASS</h3>
                    <div class="path-label">Public datasets · Web application</div>
                    <p>
                        Convert Gene Expression Omnibus studies into harmonized expression matrices and metadata prepared for downstream COMPASS analysis.
                    </p>
                    <div class="mini-list" aria-label="GEO-2-COMPASS capabilities">
                        <span>GEO accession input</span>
                        <span>Probe mapping</span>
                        <span>Annotation harmonization</span>
                    </div>
                    <a class="card-action" href="https://geo2compass.precsn.com/" target="_blank" rel="noopener noreferrer">
                        Launch GEO-2-COMPASS&nbsp; ↗
                    </a>
                </article>

                <article class="path-card soon">
                    <div class="card-status coming">Coming soon</div>
                    <h3>RAW-2-COMPASS</h3>
                    <div class="path-label">Custom datasets · Standalone software</div>
                    <p>
                        Future downloadable software for preparing your own raw sequencing outputs locally and producing standardized COMPASS-ready files.
                    </p>
                    <div class="mini-list" aria-label="Planned RAW-2-COMPASS formats">
                        <span>Raw count files</span>
                        <span>FASTQ / SRA</span>
                        <span>Common RNA-seq outputs</span>
                    </div>
                    <button class="disabled-action" type="button" disabled aria-disabled="true" title="RAW-2-COMPASS is coming soon">
                        Download RAW-2-COMPASS <small>Coming soon</small>
                    </button>
                </article>
            </div>
        </section>

        <section class="workflow-wrap" id="workflow">
            <div class="workflow-top">
                <h2>From dataset to a consistent COMPASS-ready input.</h2>
                <p>COMPASS Prep makes the transformation explicit, traceable, and approachable—without requiring a custom preprocessing script.</p>
            </div>
            <div class="workflow">
                <article class="step">
                    <div class="step-number">01 · SELECT</div>
                    <h3>Choose your data source</h3>
                    <p>Start with a public GEO accession today, or use your own dataset when the standalone application becomes available.</p>
                </article>
                <article class="step">
                    <div class="step-number">02 · PREPARE</div>
                    <h3>Harmonize the inputs</h3>
                    <p>Standardize identifiers, annotations, sample context, and expression values through a defined workflow.</p>
                </article>
                <article class="step">
                    <div class="step-number">03 · CONTINUE</div>
                    <h3>Move into COMPASS</h3>
                    <p>Export a clean expression matrix and supporting metadata ready for deterministic biological analysis.</p>
                </article>
            </div>
        </section>

        <section class="section" aria-labelledby="principles-title">
            <p class="section-kicker">Built for scientific clarity</p>
            <h2 class="section-title" id="principles-title">A dependable front door to the COMPASS ecosystem.</h2>
            <div class="principles">
                <div class="principle">
                    <strong>Standardized</strong>
                    <span>Consistent outputs across supported data sources.</span>
                </div>
                <div class="principle">
                    <strong>Deterministic</strong>
                    <span>Defined transformations designed for repeatable analysis.</span>
                </div>
                <div class="principle">
                    <strong>Reproducible</strong>
                    <span>A clear path from source dataset to prepared input.</span>
                </div>
                <div class="principle">
                    <strong>Coding-free</strong>
                    <span>Accessible workflows without scripting prerequisites.</span>
                </div>
            </div>
        </section>

        <footer class="site-footer">
            <div><span class="footer-brand">COMPASS Prep™</span> · Data harmonization and preprocessing</div>
            <div class="footer-links">
                <a href="#top">Back to top</a>
                <a href="https://geo2compass.precsn.com/" target="_blank" rel="noopener noreferrer">GEO-2-COMPASS ↗</a>
            </div>
        </footer>
    </main>
    """
)
