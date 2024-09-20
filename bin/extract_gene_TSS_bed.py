#!/usr/bin/env python
"""
Extract transcription start sites (TSSs) of protein-coding genes from a GTF file
and write out a BED file with 0.5Mbp of padding around each TSS.

Usage:
    extract_gene_TSS_bed.py <gtf_input> <chrom_sizes_input> <bed_output>

Example:
    extract_gene_TSS_bed.py /data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/hg19/gencode.v19.annotation.gtf /data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/indexes/hg19.fa.sizes hg19_TSS_padded.bed
"""
import sys


def extract_gene_tss(gtf_file):
    for line in gtf_file:
        line_tab = line.strip().split("\t")
        if len(line_tab) == 9:
            (
                seqname,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                attribute,
            ) = line_tab
            if feature == "gene":
                attr_dict = {
                    attr.strip().split(" ")[0]: attr.strip().split(" ")[1].strip('"')
                    for attr in attribute.strip(";").split(";")
                }
                if attr_dict["gene_type"] == "protein_coding":
                    tss = None
                    if strand == "+":
                        tss = int(start)
                    elif strand == "-":
                        tss = int(end)

                    yield seqname, tss, strand, attr_dict["gene_name"]


def main(gtf_filename, chrom_sizes_filename, bed_filename, padding=500000):
    padding = int(padding)
    with open(chrom_sizes_filename, "r") as infile:
        chrom_sizes = {line.split()[0]: int(line.split()[1]) for line in infile}
    with open(gtf_filename, "r") as gtf_file:
        with open(bed_filename, "w") as bed_file:
            for seqname, tss, strand, gene_name in extract_gene_tss(gtf_file):
                # make sure padding doesn't cause TSS to go over ends of chromosome
                start = tss - padding if tss - padding > 0 else 1
                end = (
                    tss + padding
                    if tss + padding <= chrom_sizes[seqname]
                    else chrom_sizes[seqname]
                )
                # write TSS ± padding to BED format
                line = [seqname, str(start), str(end), gene_name, ".", strand]
                bed_file.write("\t".join(line) + "\n")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise ValueError(
            "3 positional arguments are required: GTF filename for input, chromosome sizes filename for input, and BED filename for output"
        )
    main(*sys.argv[1:])
