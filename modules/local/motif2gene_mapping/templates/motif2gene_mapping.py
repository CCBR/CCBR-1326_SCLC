#!/usr/bin/env python
"""
Map Jaspar TF motifs to the genes that encode them
"""

import sys


def main(gtf_file, motifs_file, output_file):
    quote = '"'
    tab = "\\t"
    newline = "\\n"
    protein2gene = dict()
    with open(gtf_file, "r") as infile:
        for l in infile:
            l = l.strip().split()
            if len(l) < 3 or l[2] != "gene":
                continue
            gid = ""
            gname = ""
            for i, x in enumerate(l):
                #    print("##"+x+"##",i)
                if x == "gene_id":
                    gid = i + 1
                if x == "gene_name":
                    gname = i + 1
            if gid != "" and gname != "":
                x = l[gid].strip(";").split(".")[0] + quote
                x = x.upper()
                y = l[gname].strip(";")
                if not y in protein2gene:
                    protein2gene[y] = x
    with open(output_file, "w") as outfile:
        with open(motifs_file, "r") as infile:
            for l in infile:
                if not l.startswith(">"):
                    continue
                motif = l.strip().split()[-1]
                umotif = quote + motif.upper() + quote
                if not umotif in protein2gene:
                    continue
                outfile.write(
                    tab.join(
                        list(
                            map(lambda x: x.strip(quote), [motif, protein2gene[umotif]])
                        )
                    )
                    + newline
                )


if __name__ == "__main__":
    main("${gtf}", "${pfm}", "${output_txt}")
