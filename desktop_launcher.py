"""Launcher for the GeneWalk Desktop app.

PyInstaller bundles this script as the executable entry point.  When the user
double-clicks the resulting .exe / .app, this script starts a local Streamlit
server and opens the browser to the desktop app interface.

For development (no PyInstaller), just run:  streamlit run desktop.py
"""

import os
import sys
import threading
import time
import webbrowser


def _open_browser(port: int, delay: float = 3.0):
    """Open the default browser after a short delay to let the server start."""
    time.sleep(delay)
    webbrowser.open(f"http://localhost:{port}")


def main():
    port = 8501

    # When running as a frozen PyInstaller bundle, sys._MEIPASS points to
    # the temporary directory where PyInstaller extracted the bundled files.
    if getattr(sys, "frozen", False):
        bundle_dir = sys._MEIPASS
        app_path = os.path.join(bundle_dir, "desktop.py")

        # Point Streamlit at the bundled config
        config_dir = os.path.join(bundle_dir, ".streamlit")
        if os.path.isdir(config_dir):
            os.environ["STREAMLIT_CONFIG_DIR"] = config_dir

        # Ensure Streamlit can find its static frontend assets
        static_dir = os.path.join(bundle_dir, "streamlit", "static")
        if os.path.isdir(static_dir):
            os.environ["STREAMLIT_STATIC_DIR"] = static_dir
    else:
        app_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "desktop.py")

    sys.argv = [
        "streamlit", "run", app_path,
        f"--server.port={port}",
        "--server.headless=false",
        "--browser.gatherUsageStats=false",
        "--global.developmentMode=false",
    ]

    # Open browser in a background thread (fallback in case Streamlit's
    # own browser-open doesn't work in the frozen environment)
    threading.Thread(target=_open_browser, args=(port,), daemon=True).start()

    from streamlit.web.cli import main as st_main
    st_main()


if __name__ == "__main__":
    main()
