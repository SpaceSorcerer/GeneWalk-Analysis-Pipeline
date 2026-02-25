"""Custom CSS styles for the GeneWalk Analysis Pipeline app."""


def get_custom_css() -> str:
    """Return the full custom CSS stylesheet as a string."""
    return """
<style>
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
    background: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 10px;
    padding: 1rem 1.25rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04);
    transition: box-shadow 0.2s, transform 0.2s;
}
div[data-testid="stMetric"]:hover {
    box-shadow: 0 4px 12px rgba(0,0,0,0.08);
    transform: translateY(-1px);
}
div[data-testid="stMetric"] label {
    font-size: 0.78rem;
    font-weight: 600;
    color: #64748b;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}
div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
    font-size: 1.6rem;
    font-weight: 700;
    color: #1a365d;
}

/* ------------------------------------------------------------------ */
/* Step cards on landing page                                          */
/* ------------------------------------------------------------------ */
.step-card {
    background: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 12px;
    padding: 1.75rem;
    height: 100%;
    transition: box-shadow 0.2s, transform 0.2s;
    position: relative;
    overflow: hidden;
}
.step-card:hover {
    box-shadow: 0 8px 24px rgba(0,0,0,0.08);
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
    color: #1e293b;
    margin: 0 0 0.5rem 0;
}
.step-card p {
    font-size: 0.9rem;
    color: #64748b;
    line-height: 1.6;
    margin: 0;
}

/* ------------------------------------------------------------------ */
/* Tabs                                                                */
/* ------------------------------------------------------------------ */
.stTabs [data-baseweb="tab-list"] {
    gap: 0;
    background: #f8fafc;
    border-radius: 10px;
    padding: 4px;
    border: 1px solid #e2e8f0;
}
.stTabs [data-baseweb="tab"] {
    border-radius: 8px;
    padding: 0.6rem 1.2rem;
    font-size: 0.85rem;
    font-weight: 500;
    color: #64748b;
    background: transparent;
    border: none;
    transition: all 0.15s;
}
.stTabs [aria-selected="true"] {
    background: #ffffff !important;
    color: #1a365d !important;
    font-weight: 600 !important;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
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
section[data-testid="stSidebar"] {
    background: #f8fafc;
    border-right: 1px solid #e2e8f0;
}
section[data-testid="stSidebar"] .stMarkdown h1 {
    font-size: 1.3rem;
    font-weight: 700;
    color: #1a365d;
}
section[data-testid="stSidebar"] .stMarkdown h2,
section[data-testid="stSidebar"] .stMarkdown h3 {
    font-size: 0.85rem;
    font-weight: 600;
    color: #475569;
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
    background: #f8fafc;
    border: 1px solid #e2e8f0;
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
    border: 1px solid #e2e8f0;
}

/* ------------------------------------------------------------------ */
/* Download button                                                     */
/* ------------------------------------------------------------------ */
.stDownloadButton > button {
    border-radius: 8px;
    font-weight: 500;
    border: 1px solid #e2e8f0;
    transition: all 0.2s;
}
.stDownloadButton > button:hover {
    border-color: #4A90D9;
    color: #4A90D9;
}

/* ------------------------------------------------------------------ */
/* Section header                                                      */
/* ------------------------------------------------------------------ */
.section-header {
    font-size: 0.75rem;
    font-weight: 600;
    color: #94a3b8;
    text-transform: uppercase;
    letter-spacing: 1px;
    margin-bottom: 0.75rem;
    padding-bottom: 0.5rem;
    border-bottom: 2px solid #e2e8f0;
}

/* ------------------------------------------------------------------ */
/* Tooltip-style info boxes                                            */
/* ------------------------------------------------------------------ */
.info-tip {
    background: #eff6ff;
    border-left: 3px solid #4A90D9;
    border-radius: 0 8px 8px 0;
    padding: 0.75rem 1rem;
    font-size: 0.85rem;
    color: #1e40af;
    margin: 0.5rem 0;
}

/* ------------------------------------------------------------------ */
/* Footer                                                              */
/* ------------------------------------------------------------------ */
.app-footer {
    text-align: center;
    padding: 1.5rem 0 0.5rem 0;
    color: #94a3b8;
    font-size: 0.8rem;
    border-top: 1px solid #e2e8f0;
    margin-top: 2rem;
}
.app-footer a {
    color: #4A90D9;
    text-decoration: none;
}
.app-footer a:hover {
    text-decoration: underline;
}
</style>
"""
