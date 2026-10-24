"""101 extract a gene-panel TPM table from AoU blood RNA-seq (VM side; PRIVATE per-person output stays on the VM).

Streams the RNA-SeQC 2 gene_tpm GCT file(s) (rows genes, columns samples) in chunks and keeps a panel: KIR genes, HLA and MHC genes (classical, E/F/G, MICA/B),
C4A/C4B, complement, NK / T / B / myeloid lineage markers (composition proxies) and an interferon-stimulated-gene set. Maps sample -> research_id with the
RNA manifest. Output: <out>/rna_panel_tpm.private.tsv (person x panel gene TPM; mode 600). Downstream tests (102) join it with long-read features and report
only aggregates. Gene names are matched to the GCT 'Description' column.
"""
import argparse
import csv
import glob
import gzip
import os
import sys

import pandas as pd

PANEL = (
    # KIR
    ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL4", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DL2", "KIR3DL3", "KIR3DS1", "KIR2DP1", "KIR3DP1"]
    # HLA / MHC
    + ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "MICA", "MICB", "TAP1", "TAP2", "PSMB8", "PSMB9", "B2M"]
    # complement
    + ["C4A", "C4B", "C2", "CFB", "C3", "C1QA", "C1QB", "C1QC"]
    # lineage / composition markers
    + ["NCAM1", "NKG7", "GNLY", "KLRF1", "KLRD1", "KLRC1", "KLRK1", "GZMB", "PRF1", "CD3E", "CD3D", "CD4", "CD8A", "CD8B", "CD19", "MS4A1", "CD79A", "CD14", "LYZ", "FCGR3A", "PPBP", "HBB"]
    # interferon-stimulated genes
    + ["ISG15", "IFI44L", "IFI44", "MX1", "MX2", "IFIT1", "IFIT3", "OAS1", "OAS3", "RSAD2", "STAT1"]
    # housekeeping / negative controls
    + ["ACTB", "GAPDH", "TBP", "B2M", "PPIA", "RPL13A", "HPRT1", "UBC", "GUSB", "PGK1", "RPLP0", "TFRC", "RPS18", "YWHAZ"]
    # unfolded-protein-response / ER-associated degradation / antigen-processing / IL-23 axis (HLA-B*27 misfolding hypothesis)
    + ["HSPA5", "DDIT3", "XBP1", "ATF6", "ATF4", "EIF2AK3", "ERN1", "SEL1L", "SYVN1", "HERPUD1", "EDEM1", "DNAJB9", "PDIA4", "PDIA6", "CALR", "CANX",
       "ERAP1", "ERAP2", "TAPBP", "IL23A", "IL23R", "IL17A", "TNF", "IFNG"]
    # NK receptors for the HLA-E / NKG2A tone hypothesis
    + ["KLRC2", "KLRB1", "NCR1", "NCR3", "IL7R", "CD2", "TRDC"]
)


def log(m):
    print("[101] " + m, file=sys.stderr, flush=True)


def run(a):
    out = os.path.expanduser(a.out_dir)
    os.makedirs(out, exist_ok=True)
    root = os.path.expanduser(a.rna_root)
    files = sorted(glob.glob(os.path.join(root, "rnaseqc2", "*gene_tpm*.gct.gz")))
    log("gct files: %d" % len(files))
    sample2person = {}
    with open(os.path.expanduser(a.manifest)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sample2person.setdefault(r["sampleid"], r["research_id"])
    want = set(PANEL)
    rows = {}
    for path in files:
        with gzip.open(path, "rt") as f:
            f.readline()
            f.readline()
            header = f.readline().rstrip("\n").split("\t")
        samples = header[2:]
        got = {}
        for ch in pd.read_csv(path, sep="\t", skiprows=2, chunksize=3000, compression="gzip"):
            m = ch.iloc[:, 1].astype(str).isin(want)
            for _, r in ch[m].iterrows():
                got.setdefault(str(r.iloc[1]), r.iloc[2:].to_numpy(dtype=float))   # first occurrence per gene symbol
        log("%s: %d samples, %d panel genes found" % (os.path.basename(path), len(samples), len(got)))
        for j, s in enumerate(samples):
            p = sample2person.get(s)
            if p is None:
                continue
            rows[p] = {g: v[j] for g, v in got.items()}
    genes = [g for g in dict.fromkeys(PANEL)]
    pth = os.path.join(out, "rna_panel_tpm.private.tsv")
    with open(pth, "w") as f:
        f.write("person_id\t" + "\t".join(genes) + "\n")
        for p, d in rows.items():
            f.write(p + "\t" + "\t".join("%.4f" % d[g] if g in d else "" for g in genes) + "\n")
    os.chmod(pth, 0o600)
    log("wrote %d persons x %d panel genes" % (len(rows), len(genes)))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    base = "~/workspace/vwb-aou-datasets-controlled-v9/v9/multiomics/rnaseq"
    p.add_argument("--rna-root", default=base)
    p.add_argument("--manifest", default=base + "/manifest.tsv")
    p.add_argument("--out-dir", default="~/s09")
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
