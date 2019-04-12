#!/usr/bin/env python
from deeptoolsintervals import GTF
import py2bit
import argparse

def hasGenes(chrom, start, end, genes, flank):
    """
    Returns True if there is at least 1 gene within `flank` of the specified region
    """
    o = genes.findOverlaps(chrom, max(0, start - flank), end + flank)
    return len(o) > 0


def splitByGenes(chrom, start, end, genes, flank):
    """
    Split a region into chunks not within flank of a gene
    """
    #print("original range {}:{}-{}".format(chrom, start, end))
    o = genes.findOverlaps(chrom, max(0, start - flank), end + flank)
    o = [(max(x[0] - flank, 0), x[1] + flank) for x in o]
    #print("overlaps: {}".format(o))
    if len(o) == 0:
        return [[start, end]]

    # Merge
    o.sort()
    oFinal = []
    Gstart, Gend = o[0]
    for reg in o[1:]:
        if reg[0] <= Gend:
            if reg[1] > Gend:
                Gend = reg[1]
        else:
            oFinal.append([Gstart, Gend])
            Gstart, Gend = reg
    oFinal.append([Gstart, Gend])
    #print("final overlaps {}".format(oFinal))

    out = []
    for reg in oFinal:
        if reg[0] > start:
            out.append([start, reg[0]])
        start = reg[1]
    if start < end:
        out.append([start, end])
    #print("final: {}".format(out))

    return out


def highN(chrom, start, end, tb, threshold):
    """
    Returns True if the N content in the region is above the threshold
    """
    bases = tb.bases(chrom, start, end)
    N = 1. - sum(list(bases.values()))
    return N > threshold

    
def main():
    parser = argparse.ArgumentParser(add_help=True, description="Bourbon finds contiguous regions without repeats (low peat content) of a minimum size and without genes within some distance. Output is written to the terminal. Note that this program currently ignores the ends of chromosomes.")
    parser.add_argument("rmsk", help="Repeat masker file")
    parser.add_argument("gtf", help="GTF file")
    parser.add_argument("tbit", help="2bit file")
    parser.add_argument("--minimumProof", type=int, default=15000, help="Minimum size of a repeat-free region (default %(default)s)")
    parser.add_argument("--wobble", type=int, default=5000, help="Ensure no genes are within this distance of a region of interest (default %(default)s)")
    parser.add_argument("--legalBAC", type=float, default=0.01, help="Maximum N content (default %(default)s)")
    args = parser.parse_args()

    # Produce a header
    print("Chromosome\tStart\tEnd")

    genes = GTF(args.gtf)
    rmsk = open(args.rmsk)
    tb = py2bit.open(args.tbit)

    lastChrom = None
    lastEnd = 0
    for line in rmsk:
        if line.startswith("#"):
            continue
        cols = line.strip().split()
        chrom = cols[5]
        start = int(cols[6]) - 1
        end = int(cols[7])
        if chrom == lastChrom:
            if start - lastEnd >= args.minimumProof:
                ROIstart = lastEnd
                ROIend = start
                blocks = splitByGenes(chrom, ROIstart, ROIend, genes, args.wobble)
                for block in blocks:
                    if block[1] - block[0] < args.minimumProof:
                        continue
                    if not highN(chrom, block[0], block[1], tb, args.legalBAC):
                        print("{}\t{}\t{}".format(chrom, block[0], block[1]))
        lastChrom = chrom
        lastEnd = end

    rmsk.close()
    tb.close()

if __name__ == "__main__":
    main()
