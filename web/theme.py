"""Shared Gradio theme + CSS for ChiralFold Web / Hugging Face Space.

Designed for high contrast in both OS light and dark modes. Gradio 5 follows
``prefers-color-scheme`` and a ``.dark`` class; navy header text on a dark
surface became unreadable. This stylesheet forces a light, readable chrome
and locks explicit dark-on-light text colors (no ``color: inherit`` traps).
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

# Hex tokens used by CSS and by contrast tests (keep in sync).
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

CUSTOM_CSS = f"""
/* ---- Design tokens (WCAG-friendly on light surfaces) ---- */
:root {{
  --cf-ink: {TOKENS["ink"]};
  --cf-body: {TOKENS["body"]};
  --cf-muted: {TOKENS["muted"]};
  --cf-teal: {TOKENS["teal"]};
  --cf-teal-bright: #0d9488;
  --cf-surface: {TOKENS["surface"]};
  --cf-surface-2: {TOKENS["surface_2"]};
  --cf-border: {TOKENS["border"]};
  --cf-ok: {TOKENS["ok"]};
  --cf-warn: {TOKENS["warn"]};
  --cf-bad: {TOKENS["bad"]};
  --cf-badge-bg: {TOKENS["badge_bg"]};
  --cf-badge-fg: {TOKENS["badge_fg"]};
  --cf-badge-border: {TOKENS["badge_border"]};
}}

/* Force a light, readable app even when Gradio/OS is in dark mode */
html, body, .gradio-container,
.dark, .dark .gradio-container,
.dark body, body.dark {{
  color-scheme: light !important;
}}

html, body {{
  background: {TOKENS["page_bg"]} !important;
  color: var(--cf-body) !important;
}}

.gradio-container,
.dark .gradio-container {{
  font-family: "IBM Plex Sans", "Source Sans 3", "Segoe UI", Helvetica, Arial, sans-serif !important;
  color: var(--cf-body) !important;
  background:
    radial-gradient(ellipse 70% 45% at 8% 0%, rgba(15, 118, 110, 0.12), transparent 55%),
    radial-gradient(ellipse 50% 35% at 92% 5%, rgba(14, 116, 144, 0.10), transparent 50%),
    linear-gradient(180deg, #f8fafc 0%, #eef2f7 100%) !important;
  max-width: 1100px !important;
  margin: 0 auto !important;
}}

/* Explicit text colors — never rely on inherit from Gradio dark tokens */
.gradio-container,
.gradio-container .prose,
.gradio-container .markdown,
.gradio-container label,
.gradio-container span,
.gradio-container p,
.gradio-container li,
.dark .gradio-container,
.dark .gradio-container .prose,
.dark .gradio-container .markdown,
.dark .gradio-container label,
.dark .gradio-container span,
.dark .gradio-container p,
.dark .gradio-container li {{
  color: var(--cf-body) !important;
}}

.gradio-container h1,
.gradio-container h2,
.gradio-container h3,
.gradio-container .prose h1,
.gradio-container .prose h2,
.gradio-container .prose h3,
.dark .gradio-container h1,
.dark .gradio-container h2,
.dark .gradio-container h3 {{
  color: var(--cf-ink) !important;
}}

/* Inputs / buttons stay high-contrast */
.gradio-container input,
.gradio-container textarea,
.gradio-container select,
.dark .gradio-container input,
.dark .gradio-container textarea,
.dark .gradio-container select {{
  background: var(--cf-surface) !important;
  color: var(--cf-ink) !important;
  border-color: var(--cf-border) !important;
}}

.gradio-container button,
.dark .gradio-container button {{
  color: var(--cf-ink) !important;
}}

.gradio-container button.primary,
.gradio-container .primary,
.dark .gradio-container button.primary,
.dark .gradio-container .primary {{
  background: var(--cf-teal) !important;
  color: #ffffff !important;
  border: none !important;
  font-weight: 650 !important;
}}
.gradio-container button.primary:hover,
.dark .gradio-container button.primary:hover {{
  background: #0d9488 !important;
  color: #ffffff !important;
}}

/* Header */
#cf-header {{
  text-align: center;
  padding: 1.75rem 1rem 0.85rem;
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
  background: var(--cf-surface) !important;
  color: var(--cf-body) !important;
  border-color: var(--cf-border) !important;
}}

/* KPI cards in markdown results */
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

/* File download / status blocks */
.gradio-container .file-preview,
.gradio-container .file-preview-holder,
.gradio-container .upload-container,
.gradio-container .uploaded-files,
.gradio-container .file-name,
.gradio-container .file-size,
.gradio-container .wrap,
.dark .gradio-container .file-preview,
.dark .gradio-container .file-preview-holder,
.dark .gradio-container .upload-container,
.dark .gradio-container .uploaded-files,
.dark .gradio-container .wrap {{
  background: var(--cf-surface-2) !important;
  color: var(--cf-ink) !important;
  border-color: var(--cf-border) !important;
}}

/* Uploaded filename row — force readable text on light bar (Gradio 4/5) */
.gradio-container .file-preview *,
.gradio-container .upload-container *,
.gradio-container .uploaded-files *,
.gradio-container .file-preview-holder *,
.dark .gradio-container .file-preview *,
.dark .gradio-container .upload-container *,
.dark .gradio-container .uploaded-files *,
.dark .gradio-container .file-preview-holder * {{
  color: var(--cf-ink) !important;
}}

.gradio-container .file-preview span,
.gradio-container .file-preview a,
.gradio-container .file-name,
.gradio-container .file-size,
.dark .gradio-container .file-preview span,
.dark .gradio-container .file-preview a,
.dark .gradio-container .file-name,
.dark .gradio-container .file-size {{
  background: transparent !important;
  font-weight: 600 !important;
}}

/* File upload + download widgets */
.gradio-container .block.file,
.gradio-container [data-testid="file-upload"],
.dark .gradio-container .block.file,
.dark .gradio-container [data-testid="file-upload"] {{
  background: var(--cf-surface) !important;
  border-color: var(--cf-border) !important;
}}

/* Audit / error reports — no syntax-highlight line backgrounds */
.cf-report,
.cf-error,
.gradio-container pre.cf-report,
.gradio-container pre.cf-error,
.dark .gradio-container pre.cf-report,
.dark .gradio-container pre.cf-error {{
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
.gradio-container .prose pre span,
.dark .gradio-container pre.cf-report span,
.dark .gradio-container pre.cf-error span,
.dark .gradio-container .markdown pre span,
.dark .gradio-container .prose pre span {{
  background: transparent !important;
  color: var(--cf-ink) !important;
}}

/* Neutralize Gradio/highlight.js code fences in markdown results */
.gradio-container .markdown pre,
.gradio-container .prose pre,
.gradio-container .markdown pre code,
.gradio-container .prose pre code,
.dark .gradio-container .markdown pre,
.dark .gradio-container .prose pre,
.dark .gradio-container .markdown pre code,
.dark .gradio-container .prose pre code {{
  background: var(--cf-surface-2) !important;
  color: var(--cf-ink) !important;
  border: 1px solid var(--cf-border) !important;
  border-radius: 12px;
}}

.gradio-container .markdown code:not(pre code),
.gradio-container .prose code:not(pre code),
.dark .gradio-container .markdown code:not(pre code),
.dark .gradio-container .prose code:not(pre code) {{
  background: var(--cf-surface-2) !important;
  color: var(--cf-ink) !important;
  padding: 0.1rem 0.35rem;
  border-radius: 4px;
}}

/* Footer — full-width, no clipped install command */
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
    """Named (label, fg, bg, minimum_ratio) pairs that must stay readable."""
    s, s2 = TOKENS["surface"], TOKENS["surface_2"]
    page = TOKENS["page_bg"]
    badge = TOKENS["badge_bg"]
    teal = TOKENS["teal"]
    return [
        ("header_ink_on_page", TOKENS["ink"], page, 7.0),
        ("body_on_surface", TOKENS["body"], s, 7.0),
        ("muted_on_page", TOKENS["muted"], page, 4.5),
        ("kpi_label_on_surface2", TOKENS["muted"], s2, 4.5),
        ("kpi_value_on_surface2", TOKENS["ink"], s2, 7.0),
        ("kpi_ok_on_surface2", TOKENS["ok"], s2, 4.5),
        ("kpi_warn_on_surface2", TOKENS["warn"], s2, 4.5),
        ("kpi_bad_on_surface2", TOKENS["bad"], s2, 4.5),
        ("badge_fg_on_badge_bg", TOKENS["badge_fg"], badge, 4.5),
        ("primary_btn_white_on_teal", TOKENS["button_fg"], teal, 4.5),
    ]


def make_theme() -> Any:
    """Return a Soft Gradio theme with explicit high-contrast hues."""
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
        body_text_color_dark=TOKENS["body"],
        body_background_fill=TOKENS["page_bg"],
        body_background_fill_dark=TOKENS["page_bg"],
        block_background_fill=TOKENS["surface"],
        block_background_fill_dark=TOKENS["surface"],
        block_label_text_color=TOKENS["body"],
        block_label_text_color_dark=TOKENS["body"],
        block_title_text_color=TOKENS["ink"],
        block_title_text_color_dark=TOKENS["ink"],
        button_primary_background_fill=TOKENS["teal"],
        button_primary_background_fill_dark=TOKENS["teal"],
        button_primary_text_color=TOKENS["button_fg"],
        button_primary_text_color_dark=TOKENS["button_fg"],
        border_color_primary=TOKENS["border"],
        border_color_primary_dark=TOKENS["border"],
    )


def contrast_checklist() -> dict:
    """Static checks used by unit tests to lock accessibility tokens."""
    return {
        "ink": TOKENS["ink"],
        "body": TOKENS["body"],
        "surface": TOKENS["surface"],
        "forces_light": "color-scheme: light" in CUSTOM_CSS,
        "dark_override": ".dark .gradio-container" in CUSTOM_CSS,
        "no_color_inherit": "color: inherit" not in CUSTOM_CSS,
        "primary_button_white_text": "button.primary" in CUSTOM_CSS
        and "color: #ffffff" in CUSTOM_CSS,
        "ibm_plex": "IBM Plex Sans" in CUSTOM_CSS,
        "kpi_ok": TOKENS["ok"],
        "kpi_bad": TOKENS["bad"],
    }
