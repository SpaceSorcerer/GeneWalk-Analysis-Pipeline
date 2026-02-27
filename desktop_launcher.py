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


# ---------------------------------------------------------------------------
# Ensure only ONE browser tab ever opens.  In frozen (PyInstaller) bundles
# the Streamlit headless flag can be silently ignored, causing Streamlit to
# call webbrowser.open *in addition to* our own _open_browser thread.
# Monkey-patching webbrowser.open so it can only fire once eliminates the
# duplicate-tab problem regardless of what Streamlit does internally.
# ---------------------------------------------------------------------------
_browser_opened = False
_original_wb_open = webbrowser.open


def _open_browser_once(url, new=0, autoraise=True):
    global _browser_opened
    if _browser_opened:
        return
    _browser_opened = True
    _original_wb_open(url, new=new, autoraise=autoraise)


webbrowser.open = _open_browser_once
webbrowser.open_new = lambda url: _open_browser_once(url, new=1)
webbrowser.open_new_tab = lambda url: _open_browser_once(url, new=2)


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

    # Force headless mode via environment variable — this is the most
    # reliable way to prevent Streamlit from opening its own browser tab.
    # The CLI flag alone (--server.headless) can be ignored in frozen bundles.
    os.environ["STREAMLIT_SERVER_HEADLESS"] = "true"
    os.environ["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"

    sys.argv = [
        "streamlit", "run", app_path,
        f"--server.port={port}",
        "--server.headless=true",
        "--browser.gatherUsageStats=false",
        "--global.developmentMode=false",
    ]

    # Open browser once from our own thread — this is the ONLY browser
    # open; Streamlit's own open is suppressed by headless mode above.
    threading.Thread(target=_open_browser, args=(port,), daemon=True).start()

    # Belt-and-suspenders: also patch Streamlit's own open_browser function
    # so it becomes a no-op.  This covers cases where Streamlit bypasses the
    # webbrowser module or the headless flag is not respected in frozen bundles.
    from streamlit import cli_util
    cli_util.open_browser = lambda url: None

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
