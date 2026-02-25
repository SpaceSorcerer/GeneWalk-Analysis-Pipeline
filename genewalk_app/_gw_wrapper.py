"""Wrapper to run the GeneWalk CLI with a proper HTTP User-Agent header.

GeneWalk downloads several resource files (GO OBO, GOA GAF, PathwayCommons,
etc.) using ``urllib.request.urlretrieve``, which sends Python's default
User-Agent string.  Many servers now reject this with HTTP 403 Forbidden.

This wrapper installs a custom URL opener with a realistic User-Agent header
*before* delegating to the real GeneWalk CLI entry point, so every download
made by GeneWalk goes through the patched opener.

Usage (from subprocess):
    python -m genewalk_app._gw_wrapper --project my_project --genes genes.txt ...
"""

import urllib.request


def _install_opener():
    """Install a global urllib opener that sends a proper User-Agent."""
    opener = urllib.request.build_opener()
    opener.addheaders = [
        ("User-Agent",
         "GeneWalk/1.6 (genewalk-analysis-pipeline; "
         "+https://github.com/churchmanlab/genewalk)")
    ]
    urllib.request.install_opener(opener)


if __name__ == "__main__":
    _install_opener()
    from genewalk.cli import main
    main()
