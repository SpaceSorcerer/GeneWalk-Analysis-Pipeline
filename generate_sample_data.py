"""Generate a realistic sample genewalk_results.csv for demo purposes."""

import csv
import random
import sys
from pathlib import Path

random.seed(42)

GENES = [
    "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RAF1",
    "KRAS", "NRAS", "HRAS", "ARAF", "RPS6KA1", "RPS6KA3",
    "ELK1", "FOS", "JUN", "MYC", "DUSP1", "DUSP6", "SPRY2", "NF1",
]

GO_TERMS = {
    "biological_process": [
        ("GO:0000165", "MAPK cascade"),
        ("GO:0007265", "Ras protein signal transduction"),
        ("GO:0006468", "protein phosphorylation"),
        ("GO:0043406", "positive regulation of MAP kinase activity"),
        ("GO:0045893", "positive regulation of transcription, DNA-templated"),
        ("GO:0008284", "positive regulation of cell population proliferation"),
        ("GO:0006915", "apoptotic process"),
        ("GO:0007049", "cell cycle"),
        ("GO:0043410", "positive regulation of MAPK cascade"),
        ("GO:0071902", "positive regulation of protein serine/threonine kinase activity"),
        ("GO:0046777", "protein autophosphorylation"),
        ("GO:0035556", "intracellular signal transduction"),
        ("GO:0007169", "transmembrane receptor protein tyrosine kinase signaling pathway"),
        ("GO:0000187", "activation of MAPK activity"),
        ("GO:0032872", "regulation of stress-activated MAPK cascade"),
        ("GO:0045087", "innate immune response"),
        ("GO:0001934", "positive regulation of protein phosphorylation"),
        ("GO:0010628", "positive regulation of gene expression"),
        ("GO:0043524", "negative regulation of neuron apoptotic process"),
        ("GO:0048011", "neurotrophin TRK receptor signaling pathway"),
    ],
    "molecular_function": [
        ("GO:0004674", "protein serine/threonine kinase activity"),
        ("GO:0005524", "ATP binding"),
        ("GO:0004672", "protein kinase activity"),
        ("GO:0005515", "protein binding"),
        ("GO:0005525", "GTP binding"),
        ("GO:0003700", "DNA-binding transcription factor activity"),
        ("GO:0016301", "kinase activity"),
        ("GO:0004707", "MAP kinase activity"),
        ("GO:0046872", "metal ion binding"),
        ("GO:0019899", "enzyme binding"),
        ("GO:0004709", "MAP kinase kinase kinase activity"),
        ("GO:0004708", "MAP kinase kinase activity"),
    ],
    "cellular_component": [
        ("GO:0005829", "cytosol"),
        ("GO:0005634", "nucleus"),
        ("GO:0005886", "plasma membrane"),
        ("GO:0005737", "cytoplasm"),
        ("GO:0005654", "nucleoplasm"),
        ("GO:0005794", "Golgi apparatus"),
        ("GO:0016020", "membrane"),
        ("GO:0005856", "cytoskeleton"),
    ],
}

# Some genes are more related to certain GO terms -- weight them
GENE_AFFINITIES = {
    "BRAF":    {"MAPK cascade": 0.95, "protein serine/threonine kinase activity": 0.92, "protein phosphorylation": 0.90},
    "MAP2K1":  {"MAPK cascade": 0.93, "MAP kinase kinase activity": 0.91, "activation of MAPK activity": 0.88},
    "MAP2K2":  {"MAPK cascade": 0.91, "MAP kinase kinase activity": 0.89, "protein phosphorylation": 0.85},
    "MAPK1":   {"MAPK cascade": 0.94, "MAP kinase activity": 0.93, "protein phosphorylation": 0.90},
    "MAPK3":   {"MAPK cascade": 0.93, "MAP kinase activity": 0.91, "protein phosphorylation": 0.88},
    "RAF1":    {"MAPK cascade": 0.90, "MAP kinase kinase kinase activity": 0.92, "Ras protein signal transduction": 0.87},
    "KRAS":    {"Ras protein signal transduction": 0.95, "GTP binding": 0.93, "MAPK cascade": 0.85},
    "NRAS":    {"Ras protein signal transduction": 0.93, "GTP binding": 0.91, "MAPK cascade": 0.82},
    "HRAS":    {"Ras protein signal transduction": 0.92, "GTP binding": 0.90, "MAPK cascade": 0.80},
    "ARAF":    {"MAPK cascade": 0.85, "protein serine/threonine kinase activity": 0.83},
    "FOS":     {"DNA-binding transcription factor activity": 0.91, "positive regulation of transcription, DNA-templated": 0.93},
    "JUN":     {"DNA-binding transcription factor activity": 0.90, "positive regulation of transcription, DNA-templated": 0.92},
    "MYC":     {"positive regulation of cell population proliferation": 0.94, "cell cycle": 0.91, "positive regulation of gene expression": 0.90},
    "ELK1":    {"DNA-binding transcription factor activity": 0.88, "positive regulation of transcription, DNA-templated": 0.86},
    "NF1":     {"Ras protein signal transduction": 0.89, "GTP binding": 0.60, "negative regulation of neuron apoptotic process": 0.75},
    "DUSP1":   {"MAPK cascade": 0.70, "protein phosphorylation": 0.65},
    "DUSP6":   {"MAPK cascade": 0.75, "MAP kinase activity": 0.50},
    "SPRY2":   {"MAPK cascade": 0.68, "Ras protein signal transduction": 0.65},
    "RPS6KA1": {"protein serine/threonine kinase activity": 0.87, "protein phosphorylation": 0.85},
    "RPS6KA3": {"protein serine/threonine kinase activity": 0.85, "protein phosphorylation": 0.83},
}


def generate_row(gene, go_id, go_name, domain):
    affinity = GENE_AFFINITIES.get(gene, {}).get(go_name)
    if affinity is not None:
        sim = round(affinity + random.gauss(0, 0.03), 4)
        padj = round(max(1e-15, random.lognormvariate(-8, 2) * (1 - affinity)), 6)
    else:
        sim = round(random.uniform(0.05, 0.65), 4)
        padj = round(min(1.0, random.uniform(0.01, 1.0)), 6)

    sim = min(1.0, max(0.0, sim))
    sem = round(abs(random.gauss(0.03, 0.01)), 4)
    cilow = round(max(0, sim - 1.96 * sem), 4)
    ciupp = round(min(1, sim + 1.96 * sem), 4)
    global_padj = round(min(1.0, padj * random.uniform(0.8, 3.0)), 6)

    return {
        "hgnc_symbol": gene,
        "hgnc_id": f"HGNC:{random.randint(1000, 20000)}",
        "go_name": go_name,
        "go_id": go_id,
        "go_domain": domain,
        "sim": sim,
        "sem_sim": sem,
        "cilow": cilow,
        "ciupp": ciupp,
        "global_padj": global_padj,
        "gene_padj": padj,
    }


def main():
    out = Path(__file__).parent / "sample_data" / "sample_genewalk_results.csv"
    rows = []
    for gene in GENES:
        # Each gene gets ~15-30 GO term associations
        n_terms = random.randint(15, 30)
        all_terms = []
        for domain, terms in GO_TERMS.items():
            all_terms.extend([(go_id, go_name, domain) for go_id, go_name in terms])
        chosen = random.sample(all_terms, min(n_terms, len(all_terms)))
        for go_id, go_name, domain in chosen:
            rows.append(generate_row(gene, go_id, go_name, domain))

    fieldnames = [
        "hgnc_symbol", "hgnc_id", "go_name", "go_id", "go_domain",
        "sim", "sem_sim", "cilow", "ciupp", "global_padj", "gene_padj",
    ]
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {out}")


if __name__ == "__main__":
    main()
