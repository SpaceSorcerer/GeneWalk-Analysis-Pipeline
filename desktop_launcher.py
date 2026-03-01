"""Launcher for the GeneWalk Desktop app.

PyInstaller bundles this script as the executable entry point.  When the user
double-clicks the resulting .exe / .app, this script starts a local Streamlit
server and opens the browser to the desktop app interface.

For development (no PyInstaller), just run:  streamlit run desktop.py
"""

import os
import socket
import sys
import tempfile
import threading
import time
import webbrowser
from pathlib import Path

# Flag file used to detect whether a browser client has already connected
# (e.g. from restored browser tabs).  desktop.py writes this file when it
# starts rendering for a client.
_CLIENT_FLAG = Path(tempfile.gettempdir()) / ".genewalk_client_connected"


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


def _open_browser(port: int):
    """Open the browser only if no client has already connected.

    Waits for the Streamlit server to become healthy, then gives restored
    browser tabs a couple of seconds to reconnect.  If desktop.py has
    already written the client-connected flag (meaning a restored tab
    reconnected), we skip opening a new tab entirely.
    """
    import urllib.request

    # Wait for server to be healthy (up to 30 seconds)
    url = f"http://localhost:{port}/_stcore/health"
    for _ in range(60):
        try:
            urllib.request.urlopen(url, timeout=1)
            break
        except Exception:
            time.sleep(0.5)

    # Give restored browser tabs a moment to reconnect and trigger
    # desktop.py, which writes the client-connected flag file.
    time.sleep(2.0)

    if _CLIENT_FLAG.exists():
        # A client already connected (restored browser tab) — skip.
        return

    webbrowser.open(f"http://localhost:{port}")


def _resolve_app_path(app_filename: str) -> str:
    """Return the absolute path to a Streamlit app file.

    In a frozen PyInstaller bundle, the file lives inside sys._MEIPASS.
    Otherwise it is next to this launcher script.
    """
    if getattr(sys, "frozen", False):
        return os.path.join(sys._MEIPASS, app_filename)
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), app_filename)


def main():
    # Required for PyInstaller on Windows when any code uses multiprocessing.
    import multiprocessing
    multiprocessing.freeze_support()

    # Determine which app to launch.  --comparison flag launches the
    # GSEA + GeneWalk comparison pipeline; default is the standard desktop app.
    if "--comparison" in sys.argv:
        sys.argv.remove("--comparison")
        app_filename = "comparison_app.py"
        label = "GeneWalk + GSEA Comparison Pipeline"
    else:
        app_filename = "desktop.py"
        label = "GeneWalk"

    port = _find_free_port()

    # Delete the client-connected flag from any previous run so we get a
    # clean signal for this session.
    _CLIENT_FLAG.unlink(missing_ok=True)

    # When running as a frozen PyInstaller bundle, sys._MEIPASS points to
    # the temporary directory where PyInstaller extracted the bundled files.
    if getattr(sys, "frozen", False):
        # Point Streamlit at the bundled config
        bundle_dir = sys._MEIPASS
        config_dir = os.path.join(bundle_dir, ".streamlit")
        if os.path.isdir(config_dir):
            os.environ["STREAMLIT_CONFIG_DIR"] = config_dir

        # Ensure Streamlit can find its static frontend assets
        static_dir = os.path.join(bundle_dir, "streamlit", "static")
        if os.path.isdir(static_dir):
            os.environ["STREAMLIT_STATIC_DIR"] = static_dir

    app_path = _resolve_app_path(app_filename)

    print(f"Starting {label} on http://localhost:{port}")

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

    # Open browser in a background thread — but only if no restored tab
    # has already connected (see _open_browser logic above).
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

    # PYTHONUTF8 only takes effect when set *before* the Python interpreter
    # starts.  In a frozen PyInstaller bundle the interpreter is already
    # running, so the env var alone is not enough.  Monkey-patch builtins.open
    # so that text-mode file operations default to UTF-8 instead of the
    # system locale (cp1252 on Windows), which crashes GeneWalk's HTML report
    # generation on characters outside the cp1252 range.
    import builtins
    _builtin_open = builtins.open

    def _utf8_open(*args, **kwargs):
        # Leave binary-mode opens untouched.
        mode = args[1] if len(args) > 1 else kwargs.get("mode", "r")
        if "b" in mode:
            return _builtin_open(*args, **kwargs)
        # Default encoding to UTF-8 for text-mode opens when not specified.
        if "encoding" not in kwargs and len(args) <= 3:
            kwargs["encoding"] = "utf-8"
        return _builtin_open(*args, **kwargs)

    builtins.open = _utf8_open

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
