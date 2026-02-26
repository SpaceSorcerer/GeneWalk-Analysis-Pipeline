"""Build the GeneWalk Desktop executable with PyInstaller.

Prerequisites:
    pip install -r requirements-desktop.txt
    pip install pyinstaller

Usage:
    python build_desktop.py   # Build as a directory bundle (recommended)

The executable will be created in the ``dist/GeneWalk/`` directory.
"""

import os
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).parent
LAUNCHER = ROOT / "desktop_launcher.py"
DESKTOP_APP = ROOT / "desktop.py"
GENEWALK_APP_PKG = ROOT / "genewalk_app"
SAMPLE_DATA = ROOT / "sample_data"
STREAMLIT_CONFIG = ROOT / ".streamlit"

APP_NAME = "GeneWalk"


def find_package_dir(package_name: str) -> Path:
    """Locate an installed Python package's directory."""
    mod = __import__(package_name)
    return Path(mod.__file__).parent


def build():
    """Run PyInstaller to build the desktop executable."""
    streamlit_dir = find_package_dir("streamlit")
    plotly_dir = find_package_dir("plotly")

    # Streamlit needs its static frontend assets and runtime files
    streamlit_static = streamlit_dir / "static"

    # Use --onedir: Streamlit does not work reliably with --onefile because
    # it needs to locate its static assets and config files at runtime.
    cmd = [
        sys.executable, "-m", "PyInstaller",
        str(LAUNCHER),
        "--name", APP_NAME,
        "--onedir",
        "--clean",
        "--noconfirm",

        # ---- Data files ----
        "--add-data", f"{DESKTOP_APP}{os.pathsep}.",
        "--add-data", f"{GENEWALK_APP_PKG}{os.pathsep}genewalk_app",
        "--add-data", f"{SAMPLE_DATA}{os.pathsep}sample_data",
        "--add-data", f"{STREAMLIT_CONFIG}{os.pathsep}.streamlit",
        "--add-data", f"{streamlit_static}{os.pathsep}streamlit/static",

        # ---- Collect entire packages ----
        # Streamlit, Plotly, and their templates/assets must be fully
        # collected or the app will crash with missing-module errors.
        "--collect-all", "streamlit",
        "--collect-all", "plotly",
        "--collect-all", "altair",
        "--collect-all", "pydeck",
        "--collect-all", "genewalk",

        # ---- Hidden imports ----
        "--hidden-import", "streamlit",
        "--hidden-import", "streamlit.runtime.scriptrunner",
        "--hidden-import", "streamlit.runtime.scriptrunner.script_runner",
        "--hidden-import", "streamlit.web.cli",
        "--hidden-import", "streamlit.web.server",
        "--hidden-import", "plotly",
        "--hidden-import", "plotly.express",
        "--hidden-import", "plotly.graph_objects",
        "--hidden-import", "pandas",
        "--hidden-import", "pandas._libs.tslibs.timedeltas",
        "--hidden-import", "networkx",
        "--hidden-import", "numpy",
        "--hidden-import", "PIL",
        "--hidden-import", "pkg_resources",
        "--hidden-import", "pyarrow",
        "--hidden-import", "genewalk",
        "--hidden-import", "genewalk.cli",
        "--hidden-import", "genewalk.gene_lists",
        "--hidden-import", "genewalk.resources",
        # _gw_wrapper is never imported at the top level, but the
        # --run-genewalk dispatcher in desktop_launcher.py imports from it
        # when the exe is re-invoked as a GeneWalk subprocess.
        "--hidden-import", "genewalk_app._gw_wrapper",
    ]

    print(f"Building {APP_NAME}...")
    print(f"  Streamlit: {streamlit_dir}")
    print(f"  Plotly:    {plotly_dir}")
    print()

    result = subprocess.run(cmd, cwd=str(ROOT))
    if result.returncode != 0:
        print(f"\nBuild failed with return code {result.returncode}")
        sys.exit(result.returncode)

    dist_dir = ROOT / "dist" / APP_NAME
    exe_name = f"{APP_NAME}.exe" if sys.platform == "win32" else APP_NAME
    exe_path = dist_dir / exe_name

    print(f"\nBuild complete!")
    print(f"  Directory: {dist_dir}")
    print(f"  Executable: {exe_path}")
    print(
        "\nTo run: double-click the executable, or from a terminal:\n"
        f"  {exe_path}"
    )


if __name__ == "__main__":
    build()
