"""Build the GeneWalk Desktop executable with PyInstaller.

Prerequisites:
    pip install -r requirements-desktop.txt
    pip install pyinstaller

Usage:
    python build_desktop.py          # Build for your current platform
    python build_desktop.py --onedir # Build as a directory (faster, easier to debug)

The executable will be created in the ``dist/`` directory.
"""

import os
import shutil
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


def find_streamlit_dir() -> Path:
    """Locate the installed Streamlit package directory."""
    import streamlit
    return Path(streamlit.__file__).parent


def build(onedir: bool = False):
    """Run PyInstaller to build the desktop executable."""
    streamlit_dir = find_streamlit_dir()

    # Streamlit needs its static assets at runtime
    streamlit_static = streamlit_dir / "static"

    cmd = [
        sys.executable, "-m", "PyInstaller",
        str(LAUNCHER),
        "--name", APP_NAME,
        "--clean",
        "--noconfirm",

        # ---- Data files ----
        # The desktop Streamlit app
        "--add-data", f"{DESKTOP_APP}{os.pathsep}.",
        # The genewalk_app package
        "--add-data", f"{GENEWALK_APP_PKG}{os.pathsep}genewalk_app",
        # Sample data
        "--add-data", f"{SAMPLE_DATA}{os.pathsep}sample_data",
        # Streamlit config
        "--add-data", f"{STREAMLIT_CONFIG}{os.pathsep}.streamlit",
        # Streamlit's static assets (HTML, JS, CSS for the web UI)
        "--add-data", f"{streamlit_static}{os.pathsep}streamlit/static",

        # ---- Hidden imports ----
        # Streamlit and its dependencies that PyInstaller misses
        "--hidden-import", "streamlit",
        "--hidden-import", "streamlit.runtime.scriptrunner",
        "--hidden-import", "streamlit.web.cli",
        "--hidden-import", "plotly",
        "--hidden-import", "plotly.express",
        "--hidden-import", "plotly.graph_objects",
        "--hidden-import", "pandas",
        "--hidden-import", "networkx",
        "--hidden-import", "numpy",
        "--hidden-import", "genewalk",
        "--hidden-import", "genewalk.cli",
        "--hidden-import", "genewalk.gene_lists",
        "--hidden-import", "genewalk.resources",
        "--hidden-import", "PIL",
        "--hidden-import", "pkg_resources",

        # ---- Collect submodules ----
        "--collect-submodules", "streamlit",
        "--collect-submodules", "plotly",
        "--collect-submodules", "genewalk",
    ]

    if onedir:
        cmd.append("--onedir")
    else:
        cmd.append("--onefile")

    print(f"Building {APP_NAME}...")
    print(f"  Mode: {'onedir' if onedir else 'onefile'}")
    print(f"  Streamlit: {streamlit_dir}")
    print()

    result = subprocess.run(cmd, cwd=str(ROOT))
    if result.returncode != 0:
        print(f"\nBuild failed with return code {result.returncode}")
        sys.exit(result.returncode)

    dist_dir = ROOT / "dist"
    if onedir:
        exe_path = dist_dir / APP_NAME
        print(f"\nBuild complete!  Directory: {exe_path}")
    else:
        # Name varies by platform
        exe_name = f"{APP_NAME}.exe" if sys.platform == "win32" else APP_NAME
        exe_path = dist_dir / exe_name
        print(f"\nBuild complete!  Executable: {exe_path}")

    print(
        "\nTo run the app, double-click the executable or run it from the "
        "command line."
    )


if __name__ == "__main__":
    onedir = "--onedir" in sys.argv
    build(onedir=onedir)
