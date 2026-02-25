"""Wrapper to run the GeneWalk CLI with a proper HTTP User-Agent header.

GeneWalk downloads several resource files (GO OBO, GOA GAF, PathwayCommons,
etc.) using ``urllib.request.urlretrieve``, which sends Python's default
User-Agent string.  Many servers now reject this with HTTP 403 Forbidden.

This wrapper installs a custom URL opener with a realistic User-Agent header
*before* delegating to the real GeneWalk CLI entry point, so every download
made by GeneWalk goes through the patched opener.

It also patches the GeneMapper class to handle changes in the MGI_EntrezGene.rpt
file format from JAX, where some rows now have fewer columns than expected.

Usage (from subprocess):
    python -m genewalk_app._gw_wrapper --project my_project --genes genes.txt ...
"""

import csv
import logging
import os
import urllib.request

# Force UTF-8 for all file I/O. On Windows, Python defaults to cp1252 which
# cannot encode certain Unicode characters used in GeneWalk's HTML reports.
os.environ["PYTHONUTF8"] = "1"

logger = logging.getLogger("genewalk.gene_lists")


def _install_opener():
    """Install a global urllib opener that sends a proper User-Agent."""
    opener = urllib.request.build_opener()
    opener.addheaders = [
        ("User-Agent",
         "GeneWalk/1.6 (genewalk-analysis-pipeline; "
         "+https://github.com/churchmanlab/genewalk)")
    ]
    urllib.request.install_opener(opener)


def _patch_gene_mapper():
    """Patch GeneMapper.__init__ to handle MGI_EntrezGene.rpt format changes.

    The JAX MGI_EntrezGene.rpt file has changed over time and some rows now
    contain fewer than 9 tab-separated columns.  The upstream GeneWalk code
    unconditionally accesses ``row[8]`` (the Entrez ID column), which raises
    ``IndexError`` on short rows.  This patch skips rows that don't have
    enough columns instead of crashing.
    """
    import genewalk.gene_lists as gl

    _original_init = gl.GeneMapper.__init__

    def _patched_init(self, resource_manager):
        self.resource_manager = resource_manager
        self.hgnc_file = self.resource_manager.get_hgnc()
        self.mgi_entrez_file = self.resource_manager.get_mgi_entrez()

        # Process the MGI-Entrez mapping file with robust row handling
        self.entrez_to_mgi = {}
        with open(self.mgi_entrez_file, "r") as fh:
            csvreader = csv.reader(fh, delimiter="\t")
            for row in csvreader:
                if len(row) < 9:
                    continue
                mgi_raw = row[0]
                entrez = row[8].strip()
                if not entrez:
                    continue
                # Remove "MGI:" prefix
                mgi = mgi_raw[4:] if mgi_raw.startswith("MGI:") else mgi_raw
                self.entrez_to_mgi[entrez] = mgi

        # Delegate the rest of __init__ (HGNC processing) to the original.
        # We need to re-run the HGNC section.  The simplest safe approach is
        # to inline the remaining initialization that the original performs
        # after the MGI block.
        import re

        self.hgnc_id_to_name = {}
        self.hgnc_name_to_id = {}
        self.hgnc_withdrawn_to_new = {}
        self.hgnc_to_uniprot = {}
        self.mgi_to_hgnc = {}
        self.rgd_to_hgnc = {}
        self.entrez_to_hgnc = {}
        self.ensembl_to_hgnc = {}
        self.prev_sym_map = {}

        with open(self.hgnc_file, "r", encoding="utf-8") as fh:
            csvreader = csv.reader(fh, delimiter="\t")
            next(csvreader)  # skip header
            for row in csvreader:
                if len(row) < 10:
                    continue
                (hgnc_id, hgnc_name, description, prev_sym_entry,
                 hgnc_status, entrez_id, uniprot_id, mgi_id, rgd_id,
                 ensembl_id) = row
                hgnc_id = hgnc_id[5:]
                if hgnc_status in {"Approved", "Entry Withdrawn"}:
                    self.hgnc_id_to_name[hgnc_id] = hgnc_name
                    self.hgnc_name_to_id[hgnc_name] = hgnc_id
                elif hgnc_status == "Symbol Withdrawn":
                    m = re.match(
                        r"symbol withdrawn, see \[HGNC:(?: ?)(\d+)\]",
                        description,
                    )
                    if m:
                        new_id = m.groups()[0]
                        self.hgnc_withdrawn_to_new[hgnc_id] = new_id
                if uniprot_id:
                    self.hgnc_to_uniprot[hgnc_id] = uniprot_id
                if entrez_id:
                    self.entrez_to_hgnc[entrez_id] = hgnc_id
                if mgi_id:
                    mgi_ids = mgi_id.split(", ")
                    for mid in mgi_ids:
                        if mid.startswith("MGI:"):
                            mid = mid[4:]
                        self.mgi_to_hgnc[mid] = hgnc_id
                if rgd_id:
                    rgd_ids = rgd_id.split(", ")
                    for rid in rgd_ids:
                        if rid.startswith("RGD:"):
                            rid = rid[4:]
                        self.rgd_to_hgnc[rid] = hgnc_id
                if prev_sym_entry:
                    prev_syms = prev_sym_entry.split(", ")
                    for prev_sym in prev_syms:
                        if prev_sym in self.prev_sym_map:
                            if isinstance(self.prev_sym_map[prev_sym], list):
                                self.prev_sym_map[prev_sym].append(hgnc_id)
                            else:
                                self.prev_sym_map[prev_sym] = [
                                    self.prev_sym_map[prev_sym],
                                    hgnc_id,
                                ]
                        else:
                            self.prev_sym_map[prev_sym] = hgnc_id
                if ensembl_id:
                    self.ensembl_to_hgnc[ensembl_id] = hgnc_id
            for old_id, new_id in self.hgnc_withdrawn_to_new.items():
                if new_id in self.hgnc_id_to_name:
                    self.hgnc_id_to_name[old_id] = self.hgnc_id_to_name[new_id]

    gl.GeneMapper.__init__ = _patched_init


if __name__ == "__main__":
    _install_opener()
    _patch_gene_mapper()
    from genewalk.cli import main
    main()
