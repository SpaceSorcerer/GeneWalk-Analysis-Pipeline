"""Custom CSS styles for the GeneWalk Analysis Pipeline app."""


def get_custom_css() -> str:
    """Return the full custom CSS stylesheet as a string."""
    return """
<style>
/* ------------------------------------------------------------------ */
/* Theme-aware CSS custom properties                                   */
/* ------------------------------------------------------------------ */
:root {
    --gw-card-bg: #ffffff;
    --gw-card-border: #e2e8f0;
    --gw-card-hover-shadow: rgba(0,0,0,0.08);
    --gw-card-shadow: rgba(0,0,0,0.04);
    --gw-text-primary: #1e293b;
    --gw-text-secondary: #64748b;
    --gw-text-muted: #94a3b8;
    --gw-text-heading: #1a365d;
    --gw-text-subheading: #475569;
    --gw-bg-secondary: #f8fafc;
    --gw-bg-tertiary: #f1f5f9;
    --gw-border-color: #e2e8f0;
    --gw-info-bg: #eff6ff;
    --gw-info-text: #1e40af;
    --gw-badge-bg: #dcfce7;
    --gw-badge-text: #166534;
    --gw-code-bg: #f1f5f9;
    --gw-code-text: #1a365d;
    --gw-tab-active-bg: #ffffff;
    --gw-tab-active-text: #1a365d;
    --gw-footer-text: #94a3b8;
    --gw-link-color: #4A90D9;
}

@media (prefers-color-scheme: dark) {
    :root {
        --gw-card-bg: #1e293b;
        --gw-card-border: #334155;
        --gw-card-hover-shadow: rgba(0,0,0,0.3);
        --gw-card-shadow: rgba(0,0,0,0.2);
        --gw-text-primary: #e2e8f0;
        --gw-text-secondary: #94a3b8;
        --gw-text-muted: #64748b;
        --gw-text-heading: #93c5fd;
        --gw-text-subheading: #cbd5e1;
        --gw-bg-secondary: #0f172a;
        --gw-bg-tertiary: #1e293b;
        --gw-border-color: #334155;
        --gw-info-bg: #172554;
        --gw-info-text: #93c5fd;
        --gw-badge-bg: #14532d;
        --gw-badge-text: #86efac;
        --gw-code-bg: #1e293b;
        --gw-code-text: #93c5fd;
        --gw-tab-active-bg: #1e293b;
        --gw-tab-active-text: #93c5fd;
        --gw-footer-text: #64748b;
        --gw-link-color: #60a5fa;
    }
}

/* ------------------------------------------------------------------ */
/* Global overrides                                                    */
/* ------------------------------------------------------------------ */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
}

/* Hide default Streamlit hamburger menu and footer for a cleaner look.
   The header itself is left visible so the native sidebar toggle arrow
   (">" in top-left) works normally across all Streamlit versions. */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}

/* Smooth scrolling */
html { scroll-behavior: smooth; }

/* Main container breathing room */
.block-container {
    padding-top: 2rem;
    padding-bottom: 2rem;
    max-width: 1200px;
}

/* ------------------------------------------------------------------ */
/* Hero / header banner                                                */
/* ------------------------------------------------------------------ */
.hero-banner {
    background: linear-gradient(135deg, #1a365d 0%, #2a6cb6 60%, #4A90D9 100%);
    border-radius: 12px;
    padding: 2.5rem 2.5rem 2rem 2.5rem;
    margin-bottom: 1.5rem;
    color: white;
    position: relative;
    overflow: hidden;
}
.hero-banner::before {
    content: "";
    position: absolute;
    top: -40%;
    right: -10%;
    width: 400px;
    height: 400px;
    background: radial-gradient(circle, rgba(255,255,255,0.08) 0%, transparent 70%);
    border-radius: 50%;
}
.hero-banner h1 {
    font-size: 2rem;
    font-weight: 700;
    margin: 0 0 0.4rem 0;
    letter-spacing: -0.5px;
}
.hero-banner p {
    font-size: 1.05rem;
    opacity: 0.9;
    margin: 0;
    font-weight: 300;
}
.hero-banner .hero-badge {
    display: inline-block;
    background: rgba(255,255,255,0.15);
    border: 1px solid rgba(255,255,255,0.25);
    border-radius: 20px;
    padding: 0.25rem 0.85rem;
    font-size: 0.75rem;
    font-weight: 500;
    margin-bottom: 0.75rem;
    letter-spacing: 0.5px;
    text-transform: uppercase;
}

/* ------------------------------------------------------------------ */
/* Metric cards                                                        */
/* ------------------------------------------------------------------ */
div[data-testid="stMetric"] {
    background: var(--gw-card-bg);
    border: 1px solid var(--gw-card-border);
    border-radius: 10px;
    padding: 1rem 1.25rem;
    box-shadow: 0 1px 3px var(--gw-card-shadow);
    transition: box-shadow 0.2s, transform 0.2s;
}
div[data-testid="stMetric"]:hover {
    box-shadow: 0 4px 12px var(--gw-card-hover-shadow);
    transform: translateY(-1px);
}
div[data-testid="stMetric"] label {
    font-size: 0.78rem;
    font-weight: 600;
    color: var(--gw-text-secondary);
    text-transform: uppercase;
    letter-spacing: 0.5px;
}
div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
    font-size: 1.6rem;
    font-weight: 700;
    color: var(--gw-text-heading);
}

/* ------------------------------------------------------------------ */
/* Step cards on landing page                                          */
/* ------------------------------------------------------------------ */
.step-card {
    background: var(--gw-card-bg);
    border: 1px solid var(--gw-card-border);
    border-radius: 12px;
    padding: 1.75rem;
    height: 100%;
    transition: box-shadow 0.2s, transform 0.2s;
    position: relative;
    overflow: hidden;
}
.step-card:hover {
    box-shadow: 0 8px 24px var(--gw-card-hover-shadow);
    transform: translateY(-2px);
}
.step-card .step-number {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 36px;
    height: 36px;
    background: linear-gradient(135deg, #4A90D9, #2a6cb6);
    color: white;
    border-radius: 50%;
    font-size: 0.9rem;
    font-weight: 700;
    margin-bottom: 0.75rem;
}
.step-card h3 {
    font-size: 1.1rem;
    font-weight: 600;
    color: var(--gw-text-primary);
    margin: 0 0 0.5rem 0;
}
.step-card p {
    font-size: 0.9rem;
    color: var(--gw-text-secondary);
    line-height: 1.6;
    margin: 0;
}

/* ------------------------------------------------------------------ */
/* Tabs                                                                */
/* ------------------------------------------------------------------ */
.stTabs [data-baseweb="tab-list"] {
    gap: 0;
    background: var(--gw-bg-secondary);
    border-radius: 10px;
    padding: 4px;
    border: 1px solid var(--gw-border-color);
}
.stTabs [data-baseweb="tab"] {
    border-radius: 8px;
    padding: 0.6rem 1.2rem;
    font-size: 0.85rem;
    font-weight: 500;
    color: var(--gw-text-secondary);
    background: transparent;
    border: none;
    transition: all 0.15s;
}
.stTabs [aria-selected="true"] {
    background: var(--gw-tab-active-bg) !important;
    color: var(--gw-tab-active-text) !important;
    font-weight: 600 !important;
    box-shadow: 0 1px 3px var(--gw-card-shadow);
}
.stTabs [data-baseweb="tab-highlight"] {
    display: none;
}
.stTabs [data-baseweb="tab-border"] {
    display: none;
}

/* ------------------------------------------------------------------ */
/* Sidebar                                                             */
/* ------------------------------------------------------------------ */
section[data-testid="stSidebar"] .stMarkdown h1 {
    font-size: 1.3rem;
    font-weight: 700;
    color: var(--gw-text-heading);
}
section[data-testid="stSidebar"] .stMarkdown h2,
section[data-testid="stSidebar"] .stMarkdown h3 {
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--gw-text-subheading);
    text-transform: uppercase;
    letter-spacing: 0.5px;
}

/* Sidebar button polish */
section[data-testid="stSidebar"] .stButton > button[kind="primary"] {
    width: 100%;
    border-radius: 8px;
    padding: 0.6rem 1rem;
    font-weight: 600;
    font-size: 0.9rem;
    transition: all 0.2s;
}

/* ------------------------------------------------------------------ */
/* Filter section                                                      */
/* ------------------------------------------------------------------ */
.filter-bar {
    background: var(--gw-bg-secondary);
    border: 1px solid var(--gw-border-color);
    border-radius: 10px;
    padding: 1.25rem;
    margin-bottom: 1rem;
}

/* ------------------------------------------------------------------ */
/* Data tables                                                         */
/* ------------------------------------------------------------------ */
.stDataFrame {
    border-radius: 10px;
    overflow: hidden;
    border: 1px solid var(--gw-border-color);
}

/* ------------------------------------------------------------------ */
/* Download button                                                     */
/* ------------------------------------------------------------------ */
.stDownloadButton > button {
    border-radius: 8px;
    font-weight: 500;
    border: 1px solid var(--gw-border-color);
    transition: all 0.2s;
}
.stDownloadButton > button:hover {
    border-color: var(--gw-link-color);
    color: var(--gw-link-color);
}

/* ------------------------------------------------------------------ */
/* Section header                                                      */
/* ------------------------------------------------------------------ */
.section-header {
    font-size: 0.75rem;
    font-weight: 600;
    color: var(--gw-text-muted);
    text-transform: uppercase;
    letter-spacing: 1px;
    margin-bottom: 0.75rem;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid var(--gw-border-color);
}

/* ------------------------------------------------------------------ */
/* Tooltip-style info boxes                                            */
/* ------------------------------------------------------------------ */
.info-tip {
    background: var(--gw-info-bg);
    border-left: 3px solid var(--gw-link-color);
    border-radius: 0 8px 8px 0;
    padding: 0.75rem 1rem;
    font-size: 0.85rem;
    color: var(--gw-info-text);
    margin: 0.5rem 0;
}

/* ------------------------------------------------------------------ */
/* Guide / tutorial section                                            */
/* ------------------------------------------------------------------ */
.guide-header {
    font-size: 1.15rem;
    font-weight: 600;
    color: var(--gw-text-primary);
    margin: 1.5rem 0 0.75rem 0;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid var(--gw-border-color);
}
.guide-section h4 {
    font-size: 0.95rem;
    font-weight: 600;
    color: var(--gw-text-heading);
    margin: 0.75rem 0 0.4rem 0;
}
.guide-section p, .guide-section li {
    font-size: 0.88rem;
    color: var(--gw-text-subheading);
    line-height: 1.65;
}
.guide-section ul {
    padding-left: 1.25rem;
    margin: 0.3rem 0;
}
.guide-section code {
    background: var(--gw-code-bg);
    padding: 0.15rem 0.4rem;
    border-radius: 4px;
    font-size: 0.82rem;
    color: var(--gw-code-text);
}
.param-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 0.85rem;
    margin: 0.5rem 0;
}
.param-table th {
    text-align: left;
    padding: 0.5rem 0.75rem;
    background: var(--gw-bg-tertiary);
    color: var(--gw-text-heading);
    font-weight: 600;
    border-bottom: 2px solid var(--gw-border-color);
}
.param-table td {
    padding: 0.5rem 0.75rem;
    border-bottom: 1px solid var(--gw-bg-tertiary);
    color: var(--gw-text-subheading);
    vertical-align: top;
}
.param-table tr:hover td {
    background: var(--gw-bg-secondary);
}
.recommend-badge {
    display: inline-block;
    background: var(--gw-badge-bg);
    color: var(--gw-badge-text);
    border-radius: 4px;
    padding: 0.1rem 0.4rem;
    font-size: 0.75rem;
    font-weight: 600;
}
.gene-example-box {
    background: var(--gw-bg-secondary);
    border: 1px solid var(--gw-border-color);
    border-radius: 8px;
    padding: 1rem 1.25rem;
    font-family: 'Courier New', monospace;
    font-size: 0.85rem;
    color: var(--gw-text-primary);
    line-height: 1.7;
    margin: 0.5rem 0;
}

/* ------------------------------------------------------------------ */
/* Footer                                                              */
/* ------------------------------------------------------------------ */
.app-footer {
    text-align: center;
    padding: 1.5rem 0 0.5rem 0;
    color: var(--gw-footer-text);
    font-size: 0.8rem;
    border-top: 1px solid var(--gw-border-color);
    margin-top: 2rem;
}
.app-footer a {
    color: var(--gw-link-color);
    text-decoration: none;
}
.app-footer a:hover {
    text-decoration: underline;
}
</style>
"""
