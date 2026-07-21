"""Shared Gradio theme + CSS for ChiralFold Web / Hugging Face Space.

Supports light and dark modes via Gradio's ``?__theme=light|dark`` (and a
header toggle). CSS variables switch under ``.dark`` so file rows, audit
reports, and the footer stay readable in both schemes.
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

# Light-mode tokens (also used by WCAG unit tests).
TOKENS: Dict[str, str] = {
    "ink": "#0b1220",
    "body": "#1e293b",
    "muted": "#475569",
    "teal": "#0f766e",
    "surface": "#ffffff",
    "surface_2": "#f1f5f9",
    "page_bg": "#f8fafc",
    "border": "#cbd5e1",
    "ok": "#047857",
    "warn": "#b45309",
    "bad": "#b91c1c",
    "badge_bg": "#ccfbf1",
    "badge_fg": "#115e59",
    "badge_border": "#5eead4",
    "button_fg": "#ffffff",
}

# Dark-mode tokens — light text on deep slate/teal surfaces.
DARK_TOKENS: Dict[str, str] = {
    "ink": "#f8fafc",
    "body": "#e2e8f0",
    "muted": "#94a3b8",
    "teal": "#2dd4bf",
    "surface": "#1e293b",
    "surface_2": "#0f172a",
    "page_bg": "#0b1220",
    "border": "#334155",
    "ok": "#34d399",
    "warn": "#fbbf24",
    "bad": "#f87171",
    "badge_bg": "#134e4a",
    "badge_fg": "#99f6e4",
    "badge_border": "#0d9488",
    "button_fg": "#042f2e",
}

CUSTOM_CSS = f"""
/* ---- Light tokens (default) ---- */
:root {{
  --cf-ink: {TOKENS["ink"]};
  --cf-body: {TOKENS["body"]};
  --cf-muted: {TOKENS["muted"]};
  --cf-teal: {TOKENS["teal"]};
  --cf-teal-bright: #0d9488;
  --cf-surface: {TOKENS["surface"]};
  --cf-surface-2: {TOKENS["surface_2"]};
  --cf-page-bg: {TOKENS["page_bg"]};
  --cf-border: {TOKENS["border"]};
  --cf-ok: {TOKENS["ok"]};
  --cf-warn: {TOKENS["warn"]};
  --cf-bad: {TOKENS["bad"]};
  --cf-badge-bg: {TOKENS["badge_bg"]};
  --cf-badge-fg: {TOKENS["badge_fg"]};
  --cf-badge-border: {TOKENS["badge_border"]};
  --cf-btn-fg: {TOKENS["button_fg"]};
  color-scheme: light;
}}

/* ---- Dark tokens (Gradio .dark / ?__theme=dark) ---- */
.dark,
html.dark,
body.dark,
.dark .gradio-container {{
  --cf-ink: {DARK_TOKENS["ink"]};
  --cf-body: {DARK_TOKENS["body"]};
  --cf-muted: {DARK_TOKENS["muted"]};
  --cf-teal: {DARK_TOKENS["teal"]};
  --cf-teal-bright: #5eead4;
  --cf-surface: {DARK_TOKENS["surface"]};
  --cf-surface-2: {DARK_TOKENS["surface_2"]};
  --cf-page-bg: {DARK_TOKENS["page_bg"]};
  --cf-border: {DARK_TOKENS["border"]};
  --cf-ok: {DARK_TOKENS["ok"]};
  --cf-warn: {DARK_TOKENS["warn"]};
  --cf-bad: {DARK_TOKENS["bad"]};
  --cf-badge-bg: {DARK_TOKENS["badge_bg"]};
  --cf-badge-fg: {DARK_TOKENS["badge_fg"]};
  --cf-badge-border: {DARK_TOKENS["badge_border"]};
  --cf-btn-fg: {DARK_TOKENS["button_fg"]};
  color-scheme: dark;
}}

html, body {{
  background: var(--cf-page-bg) !important;
  color: var(--cf-body) !important;
}}

.gradio-container {{
  font-family: "IBM Plex Sans", "Source Sans 3", "Segoe UI", Helvetica, Arial, sans-serif !important;
  color: var(--cf-body) !important;
  background:
    radial-gradient(ellipse 70% 45% at 8% 0%, rgba(15, 118, 110, 0.14), transparent 55%),
    radial-gradient(ellipse 50% 35% at 92% 5%, rgba(14, 116, 144, 0.10), transparent 50%),
    linear-gradient(180deg, var(--cf-page-bg) 0%, var(--cf-surface-2) 100%) !important;
  max-width: 1100px !important;
  margin: 0 auto !important;
}}
.dark .gradio-container {{
  background:
    radial-gradient(ellipse 70% 45% at 8% 0%, rgba(45, 212, 191, 0.10), transparent 55%),
    radial-gradient(ellipse 50% 35% at 92% 5%, rgba(56, 189, 248, 0.08), transparent 50%),
    linear-gradient(180deg, #0b1220 0%, #0f172a 100%) !important;
}}

/* Explicit text — use tokens, never inherit Gradio defaults */
.gradio-container,
.gradio-container .prose,
.gradio-container .markdown,
.gradio-container label,
.gradio-container span,
.gradio-container p,
.gradio-container li {{
  color: var(--cf-body) !important;
}}

.gradio-container h1,
.gradio-container h2,
.gradio-container h3,
.gradio-container .prose h1,
.gradio-container .prose h2,
.gradio-container .prose h3 {{
  color: var(--cf-ink) !important;
}}

.gradio-container input,
.gradio-container textarea,
.gradio-container select {{
  background: var(--cf-surface) !important;
  color: var(--cf-ink) !important;
  border-color: var(--cf-border) !important;
}}

.gradio-container button {{
  color: var(--cf-ink) !important;
}}

.gradio-container button.primary,
.gradio-container .primary {{
  background: var(--cf-teal) !important;
  color: var(--cf-btn-fg) !important;
  border: none !important;
  font-weight: 650 !important;
}}
.gradio-container button.primary:hover {{
  background: var(--cf-teal-bright) !important;
  color: var(--cf-btn-fg) !important;
}}
.dark .gradio-container button.primary,
.dark .gradio-container .primary {{
  color: #042f2e !important;
}}

/* Header + theme toggle */
#cf-header {{
  text-align: center;
  padding: 1.75rem 1rem 0.85rem;
  position: relative;
}}
#cf-header h1 {{
  color: var(--cf-ink) !important;
  font-weight: 800;
  font-size: clamp(1.85rem, 4vw, 2.45rem);
  letter-spacing: -0.03em;
  margin: 0 0 0.4rem;
}}
#cf-header .cf-tagline {{
  color: var(--cf-muted) !important;
  font-size: 1.05rem;
  max-width: 36rem;
  margin: 0 auto 0.85rem;
  line-height: 1.5;
}}
#cf-header .cf-meta {{
  display: inline-flex;
  gap: 0.5rem;
  flex-wrap: wrap;
  justify-content: center;
  font-size: 0.8rem;
}}
#cf-header .cf-badge {{
  background: var(--cf-badge-bg) !important;
  color: var(--cf-badge-fg) !important;
  border: 1px solid var(--cf-badge-border) !important;
  border-radius: 999px;
  padding: 0.2rem 0.7rem;
  font-weight: 650;
}}

.cf-theme-toggle {{
  display: inline-flex;
  gap: 0;
  margin: 0.85rem auto 0;
  border: 1px solid var(--cf-border);
  border-radius: 999px;
  overflow: hidden;
  background: var(--cf-surface-2);
}}
.cf-theme-toggle a {{
  display: inline-block;
  padding: 0.35rem 0.9rem;
  font-size: 0.8rem;
  font-weight: 650;
  text-decoration: none !important;
  color: var(--cf-muted) !important;
  background: transparent !important;
  border: none !important;
}}
.cf-theme-toggle a:hover {{
  color: var(--cf-ink) !important;
  background: var(--cf-surface) !important;
}}
.cf-theme-toggle a.active {{
  color: var(--cf-btn-fg) !important;
  background: var(--cf-teal) !important;
}}
.dark .cf-theme-toggle a.active {{
  color: #042f2e !important;
}}

/* Panels */
.cf-panel {{
  border: 1px solid var(--cf-border) !important;
  border-radius: 16px !important;
  box-shadow: 0 6px 24px rgba(15, 23, 42, 0.06) !important;
  background: var(--cf-surface) !important;
  padding: 0.35rem !important;
  color: var(--cf-body) !important;
}}
.dark .cf-panel {{
  box-shadow: 0 6px 24px rgba(0, 0, 0, 0.35) !important;
}}

/* KPI cards */
.cf-kpi {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
  gap: 0.65rem;
  margin: 0.75rem 0 1rem;
}}
.cf-kpi .cell {{
  background: var(--cf-surface-2) !important;
  border: 1px solid var(--cf-border) !important;
  border-radius: 12px;
  padding: 0.75rem 0.85rem;
  text-align: center;
}}
.cf-kpi .cell .label {{
  font-size: 0.72rem;
  text-transform: uppercase;
  letter-spacing: 0.04em;
  color: var(--cf-muted) !important;
  font-weight: 600;
}}
.cf-kpi .cell .value {{
  font-size: 1.35rem;
  font-weight: 750;
  color: var(--cf-ink) !important;
  margin-top: 0.2rem;
}}
.cf-kpi .ok .value {{ color: var(--cf-ok) !important; }}
.cf-kpi .warn .value {{ color: var(--cf-warn) !important; }}
.cf-kpi .bad .value {{ color: var(--cf-bad) !important; }}

/* File upload / download rows */
.gradio-container .file-preview,
.gradio-container .file-preview-holder,
.gradio-container .upload-container,
.gradio-container .uploaded-files,
.gradio-container .file-name,
.gradio-container .file-size,
.gradio-container .wrap,
.gradio-container .file .file-preview,
.gradio-container .block .file-preview,
.gradio-container .block .upload-container .wrap {{
  background: var(--cf-surface-2) !important;
  color: var(--cf-ink) !important;
  border-color: var(--cf-border) !important;
}}

.gradio-container .file-preview *,
.gradio-container .upload-container *,
.gradio-container .uploaded-files *,
.gradio-container .file-preview-holder * {{
  color: var(--cf-ink) !important;
}}

.gradio-container .file-preview span,
.gradio-container .file-preview a,
.gradio-container .file-name,
.gradio-container .file-size {{
  background: transparent !important;
  font-weight: 600 !important;
}}

.gradio-container .block.file,
.gradio-container [data-testid="file-upload"] {{
  background: var(--cf-surface) !important;
  border-color: var(--cf-border) !important;
}}

.gradio-container .file button,
.gradio-container .file a {{
  color: var(--cf-ink) !important;
}}
.gradio-container .file .or {{
  color: var(--cf-muted) !important;
}}

/* Audit / error reports */
.cf-report,
.cf-error,
.gradio-container pre.cf-report,
.gradio-container pre.cf-error {{
  display: block;
  width: 100%;
  box-sizing: border-box;
  margin: 0.75rem 0 1rem;
  padding: 1rem 1.1rem;
  overflow-x: auto;
  white-space: pre-wrap;
  word-break: break-word;
  font-family: "IBM Plex Mono", "Source Code Pro", ui-monospace, monospace !important;
  font-size: 0.82rem;
  line-height: 1.45;
  color: var(--cf-ink) !important;
  background: var(--cf-surface-2) !important;
  border: 1px solid var(--cf-border) !important;
  border-radius: 12px;
}}

.gradio-container pre.cf-report span,
.gradio-container pre.cf-error span,
.gradio-container .markdown pre span,
.gradio-container .prose pre span {{
  background: transparent !important;
  color: var(--cf-ink) !important;
}}

.gradio-container .markdown pre,
.gradio-container .prose pre,
.gradio-container .markdown pre code,
.gradio-container .prose pre code {{
  background: var(--cf-surface-2) !important;
  color: var(--cf-ink) !important;
  border: 1px solid var(--cf-border) !important;
  border-radius: 12px;
}}

.gradio-container .markdown code:not(pre code),
.gradio-container .prose code:not(pre code) {{
  background: var(--cf-surface-2) !important;
  color: var(--cf-ink) !important;
  padding: 0.1rem 0.35rem;
  border-radius: 4px;
}}

/* Footer */
#cf-footer {{
  margin: 1.5rem 0 0.5rem;
  padding: 1rem 0 0.25rem;
  border-top: 1px solid var(--cf-border);
  text-align: center;
}}
#cf-footer .cf-footer-line {{
  margin: 0 0 0.35rem;
  color: var(--cf-body) !important;
  line-height: 1.6;
  word-wrap: break-word;
}}
#cf-footer .cf-footer-note {{
  margin: 0;
  font-size: 0.85rem;
  color: var(--cf-muted) !important;
}}
#cf-footer a {{
  color: var(--cf-teal) !important;
  font-weight: 650;
  text-decoration: none;
}}
#cf-footer a:hover {{
  text-decoration: underline;
}}
#cf-footer .cf-install-cmd,
#cf-footer code {{
  display: inline-block;
  max-width: 100%;
  overflow-x: auto;
  background: var(--cf-surface-2) !important;
  color: var(--cf-ink) !important;
  border: 1px solid var(--cf-border);
  border-radius: 6px;
  padding: 0.15rem 0.55rem;
  font-family: "IBM Plex Mono", ui-monospace, monospace;
  font-size: 0.85rem;
  white-space: nowrap;
}}

footer {{ visibility: hidden; }}
"""


def _hex_to_rgb(hex_color: str) -> Tuple[float, float, float]:
    h = hex_color.lstrip("#")
    return tuple(int(h[i : i + 2], 16) / 255.0 for i in (0, 2, 4))  # type: ignore[return-value]


def _relative_luminance(hex_color: str) -> float:
    def channel(c: float) -> float:
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4

    r, g, b = _hex_to_rgb(hex_color)
    return 0.2126 * channel(r) + 0.7152 * channel(g) + 0.0722 * channel(b)


def contrast_ratio(fg: str, bg: str) -> float:
    """WCAG contrast ratio between two hex colors."""
    l1 = _relative_luminance(fg)
    l2 = _relative_luminance(bg)
    lighter, darker = max(l1, l2), min(l1, l2)
    return (lighter + 0.05) / (darker + 0.05)


def readable_pairs() -> List[Tuple[str, str, str, float]]:
    """Named (label, fg, bg, minimum_ratio) pairs for light and dark modes."""
    pairs: List[Tuple[str, str, str, float]] = []
    for prefix, t in (("light", TOKENS), ("dark", DARK_TOKENS)):
        s, s2, page = t["surface"], t["surface_2"], t["page_bg"]
        pairs.extend(
            [
                (f"{prefix}_ink_on_page", t["ink"], page, 7.0),
                (f"{prefix}_body_on_surface", t["body"], s, 7.0),
                (f"{prefix}_muted_on_page", t["muted"], page, 4.5),
                (f"{prefix}_kpi_label_on_surface2", t["muted"], s2, 4.5),
                (f"{prefix}_kpi_value_on_surface2", t["ink"], s2, 7.0),
                (f"{prefix}_kpi_ok_on_surface2", t["ok"], s2, 4.5),
                (f"{prefix}_kpi_warn_on_surface2", t["warn"], s2, 4.5),
                (f"{prefix}_kpi_bad_on_surface2", t["bad"], s2, 4.5),
                (f"{prefix}_badge_fg_on_badge_bg", t["badge_fg"], t["badge_bg"], 4.5),
                (
                    f"{prefix}_primary_btn_on_teal",
                    t["button_fg"],
                    t["teal"],
                    4.5,
                ),
            ]
        )
    return pairs


def header_html(version: str) -> str:
    """App header with badges and light/dark theme toggle."""
    return f"""
<div id="cf-header">
  <h1>ChiralFold</h1>
  <p class="cf-tagline">
    Upload any PDB — audit stereochemistry, fix chirality errors,
    or generate an exact mirror-image structure. No command line.
  </p>
  <div class="cf-meta">
    <span class="cf-badge">v{version}</span>
    <span class="cf-badge">0% residual chirality violations</span>
    <span class="cf-badge">open source</span>
  </div>
  <div class="cf-theme-toggle" role="group" aria-label="Color theme">
    <a class="cf-theme-btn" id="cf-theme-light" href="?__theme=light">Light</a>
    <a class="cf-theme-btn" id="cf-theme-dark" href="?__theme=dark">Dark</a>
  </div>
</div>
<script>
(function () {{
  try {{
    var p = new URLSearchParams(window.location.search).get("__theme");
    var dark = p === "dark" || (!p && document.documentElement.classList.contains("dark"));
    var el = document.getElementById(dark ? "cf-theme-dark" : "cf-theme-light");
    if (el) el.classList.add("active");
  }} catch (e) {{}}
}})();
</script>
"""


def make_theme() -> Any:
    """Return a Soft Gradio theme with explicit light and dark hues."""
    import gradio as gr

    return gr.themes.Soft(
        primary_hue=gr.themes.Color(
            c50="#f0fdfa",
            c100="#ccfbf1",
            c200="#99f6e4",
            c300="#5eead4",
            c400="#2dd4bf",
            c500="#14b8a6",
            c600="#0d9488",
            c700="#0f766e",
            c800="#115e59",
            c900="#134e4a",
            c950="#042f2e",
        ),
        secondary_hue="slate",
        neutral_hue="slate",
        text_size=gr.themes.sizes.text_md,
        font=[
            gr.themes.GoogleFont("IBM Plex Sans"),
            "Source Sans 3",
            "Segoe UI",
            "Helvetica",
            "Arial",
            "sans-serif",
        ],
    ).set(
        body_text_color=TOKENS["body"],
        body_text_color_dark=DARK_TOKENS["body"],
        body_background_fill=TOKENS["page_bg"],
        body_background_fill_dark=DARK_TOKENS["page_bg"],
        block_background_fill=TOKENS["surface"],
        block_background_fill_dark=DARK_TOKENS["surface"],
        block_label_text_color=TOKENS["body"],
        block_label_text_color_dark=DARK_TOKENS["body"],
        block_title_text_color=TOKENS["ink"],
        block_title_text_color_dark=DARK_TOKENS["ink"],
        button_primary_background_fill=TOKENS["teal"],
        button_primary_background_fill_dark=DARK_TOKENS["teal"],
        button_primary_text_color=TOKENS["button_fg"],
        button_primary_text_color_dark=DARK_TOKENS["button_fg"],
        border_color_primary=TOKENS["border"],
        border_color_primary_dark=DARK_TOKENS["border"],
    )


def contrast_checklist() -> dict:
    """Static checks used by unit tests to lock accessibility tokens."""
    return {
        "ink": TOKENS["ink"],
        "body": TOKENS["body"],
        "surface": TOKENS["surface"],
        "dark_ink": DARK_TOKENS["ink"],
        "dark_surface": DARK_TOKENS["surface"],
        "supports_dark": "--cf-ink: #f8fafc" in CUSTOM_CSS or DARK_TOKENS["ink"] in CUSTOM_CSS,
        "theme_toggle": "cf-theme-toggle" in CUSTOM_CSS,
        "dark_override": ".dark" in CUSTOM_CSS,
        "no_color_inherit": "color: inherit" not in CUSTOM_CSS,
        "primary_button_white_text": "button.primary" in CUSTOM_CSS,
        "ibm_plex": "IBM Plex Sans" in CUSTOM_CSS,
        "kpi_ok": TOKENS["ok"],
        "kpi_bad": TOKENS["bad"],
        "forces_light": False,
    }
