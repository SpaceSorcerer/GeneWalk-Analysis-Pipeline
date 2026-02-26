"""Launcher for the GeneWalk Desktop app.

PyInstaller bundles this script as the executable entry point.  When the user
double-clicks the resulting .exe / .app, this script starts a local Streamlit
server and opens the browser to the desktop app interface.

For development (no PyInstaller), just run:  streamlit run desktop.py
"""

import os
import socket
import sys
import threading
import time
import webbrowser


def _find_free_port(start: int = 8501, end: int = 8600) -> int:
    """Find the first available port in the given range."""
    for port in range(start, end):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("localhost", port))
                return port
            except OSError:
                continue
    # Fallback: let the OS pick any free port
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("localhost", 0))
        return s.getsockname()[1]


def _open_browser(port: int, delay: float = 3.0):
    """Open the default browser after a short delay to let the server start."""
    time.sleep(delay)
    webbrowser.open(f"http://localhost:{port}")


def main():
    port = _find_free_port()

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

    print(f"Starting GeneWalk on http://localhost:{port}")

    sys.argv = [
        "streamlit", "run", app_path,
        f"--server.port={port}",
        "--server.headless=true",
        "--browser.gatherUsageStats=false",
        "--global.developmentMode=false",
    ]

    # Open browser in a background thread (fallback in case Streamlit's
    # own browser-open doesn't work in the frozen environment)
    threading.Thread(target=_open_browser, args=(port,), daemon=True).start()

    from streamlit.web.cli import main as st_main
    st_main()


def _run_genewalk_subprocess():
    """Dispatch to the GeneWalk CLI when the exe is re-invoked as a subprocess.

    When ``runner.py`` spawns a GeneWalk analysis, it calls
    ``[sys.executable, "--run-genewalk", ...]``.  In the frozen bundle
    ``sys.executable`` is this very ``.exe``, so we detect the flag here
    and run GeneWalk directly instead of starting another Streamlit server.
    """
    os.environ["PYTHONUTF8"] = "1"

    # Build a clean argv for GeneWalk's argparse (strip our sentinel flag)
    sys.argv = ["genewalk"] + [a for a in sys.argv[1:] if a != "--run-genewalk"]

    # Apply the same patches that _gw_wrapper.py applies:
    # 1. Proper HTTP User-Agent so resource downloads aren't blocked
    # 2. GeneMapper patch for MGI_EntrezGene.rpt format changes
    from genewalk_app._gw_wrapper import _install_opener, _patch_gene_mapper
    _install_opener()
    _patch_gene_mapper()

    from genewalk.cli import main as gw_main
    gw_main()


if __name__ == "__main__":
    if "--run-genewalk" in sys.argv:
        _run_genewalk_subprocess()
    else:
        main()
