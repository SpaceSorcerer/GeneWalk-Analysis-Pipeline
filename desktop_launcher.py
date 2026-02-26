"""Launcher for the GeneWalk Desktop app.

PyInstaller bundles this script as the executable entry point.  When the user
double-clicks the resulting .exe / .app, this script starts a local Streamlit
server and opens the browser to the desktop app interface.

For development (no PyInstaller), just run:  streamlit run desktop.py
"""

import os
import sys


def main():
    # When running as a frozen PyInstaller bundle, ``sys._MEIPASS`` points to
    # the temporary directory where PyInstaller extracted the bundled files.
    if getattr(sys, "frozen", False):
        bundle_dir = sys._MEIPASS
        app_path = os.path.join(bundle_dir, "desktop.py")
    else:
        app_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "desktop.py")

    sys.argv = [
        "streamlit", "run", app_path,
        "--server.headless=true",
        "--browser.gatherUsageStats=false",
        "--global.developmentMode=false",
    ]

    from streamlit.web.cli import main as st_main
    st_main()


if __name__ == "__main__":
    main()
