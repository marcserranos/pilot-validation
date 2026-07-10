#!/usr/bin/env python3
"""Apply two bug fixes to ~/tools/SpecImmune, found and verified locally on 2026-07-10
while investigating Experiment C (gene-panel restriction). Neither fix is specific to
gene restriction -- both are real, pre-existing bugs in stock SpecImmune -- but the
first blocks the --HLA_fa/--HLA_exon_fa mechanism this experiment depends on.

Fix 1 (make_db.py): the --HLA_fa/--HLA_exon_fa branches reference local_release_version/
local_g_group_annotation before assigning them -- UnboundLocalError, 100% reproducible,
on every call that passes a custom reference FASTA. This is why gene-panel restriction
was "ruled out" on 2026-07-08 without actually being tried: the documented lever was
broken from the start.

Fix 2 (alignment_modules.py): cout_read_num() opens a gzipped FASTQ with plain
open(path, "rb") and counts raw compressed bytes as if they were newline-delimited
lines -- a nonsense count on any platform. Switched to gzip.open() so the count (and
the branch it drives, in subsample_fastq()) is correct. Also switched the shell-out from
`zcat` to `gzip -dc`, more portable across BSD/GNU userlands.

This patch is idempotent -- running it twice is a no-op (it checks for the fix already
being present before touching a file).

Usage: python3 patch_specimmune_for_gene_restriction.py [/path/to/SpecImmune]
Default path: ~/tools/SpecImmune
"""
import os
import sys

root = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/tools/SpecImmune")
make_db = os.path.join(root, "scripts", "make_db.py")
align_mod = os.path.join(root, "scripts", "alignment_modules.py")


def replace_once(path, old, new, label):
    with open(path) as f:
        content = f.read()
    if new in content:
        print(f"{label}: already patched, skipping")
        return
    if old not in content:
        print(f"{label}: FAILED -- expected original text not found, check manually")
        sys.exit(1)
    content = content.replace(old, new, 1)
    with open(path, "w") as f:
        f.write(content)
    print(f"{label}: patched")


replace_once(
    make_db,
    '''    else:
        local_fasta_filename = args.HLA_fa
        download_file(release_version, local_release_version)
        download_file(g_group_annotation, local_g_group_annotation)''',
    '''    else:
        local_fasta_filename = args.HLA_fa
        local_release_version = os.path.join(HLA_dir, "release_version.txt")
        local_g_group_annotation = os.path.join(HLA_dir, "hla_nom_g.txt")
        download_file(release_version, local_release_version)
        download_file(g_group_annotation, local_g_group_annotation)''',
    "make_db.py fix 1/2 (--HLA_fa branch)",
)

replace_once(
    make_db,
    '''    else:
        local_fasta_filename = args.HLA_exon_fa
        download_file(release_version, local_release_version)
        download_file(g_group_annotation, local_g_group_annotation)''',
    '''    else:
        local_fasta_filename = args.HLA_exon_fa
        local_release_version = os.path.join(HLA_exon_dir, "release_version.txt")
        local_g_group_annotation = os.path.join(HLA_exon_dir, "hla_nom_g.txt")
        download_file(release_version, local_release_version)
        download_file(g_group_annotation, local_g_group_annotation)''',
    "make_db.py fix 2/2 (--HLA_exon_fa branch)",
)

replace_once(
    align_mod,
    '''def cout_read_num(fastq):
    ## the input fastq is gziped, count the reads num in it, and return it to the variable
    f = open(fastq, "rb")
    count = 0
    for line in f:
        count += 1
    f.close()
    return count/4

def subsample_fastq(gene, read_num, my_folder):
    fastq = my_folder.reads_dir + "/" + gene + ".long_read.fq.gz"
    sub_fastq = my_folder.reads_dir + "/" + gene + ".long_read.sub.fq"
    read_count = cout_read_num(fastq)
    if read_count <= read_num:
        os.system(f"zcat {fastq} > {sub_fastq}")''',
    '''def cout_read_num(fastq):
    ## the input fastq is gziped, count the reads num in it, and return it to the variable
    import gzip
    count = 0
    with gzip.open(fastq, "rb") as f:
        for line in f:
            count += 1
    return count/4

def subsample_fastq(gene, read_num, my_folder):
    fastq = my_folder.reads_dir + "/" + gene + ".long_read.fq.gz"
    sub_fastq = my_folder.reads_dir + "/" + gene + ".long_read.sub.fq"
    read_count = cout_read_num(fastq)
    if read_count <= read_num:
        os.system(f"gzip -dc {fastq} > {sub_fastq}")''',
    "alignment_modules.py fix (gzip read-count + decompression)",
)

print("Done.")
