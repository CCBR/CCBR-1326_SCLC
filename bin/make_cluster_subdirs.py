#!/usr/bin/env python
import os


def main():
    outdir = "/data/SCLCgenomics/CCBR_analysis/ATAC_bams/clusters/"
    with open(
        "/data/SCLCgenomics/CCBR_analysis/NMF/cluster_membership.tsv", "r"
    ) as infile:
        header = next(infile)
        for line in infile:
            sample, cluster, bamfile = line.split()
            cluster_id = f"c{cluster.strip('cluster')}"
            bam_source = os.readlink(bamfile)
            os.symlink(
                bam_source,
                os.path.join(outdir, cluster_id, os.path.basename(bam_source)),
            )


if __name__ == "__main__":
    main()
