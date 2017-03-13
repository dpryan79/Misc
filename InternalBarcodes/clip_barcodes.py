#!/usr/bin/python
import sys
import os
import argparse
import gzip
try:
    import editdistance as ed
except:
    sys.exit("You must install 'editdistance' first")


def parseArgs(args=None):
    parser = argparse.ArgumentParser(description="Process RELACS reads by moving a barcode into the header and writing output to per-barcode files.")
    parser.add_argument('-b', '--buffer', help='Number of bases of "buffer" to ignore after the barcode. (default: 1)', type=int, default=1)
    parser.add_argument('sampleTable', help="""A tab-separated table with two columns: barcode	sample_name

An example is:

ACTACT	Sample1
TGACTG	Sample2
default	Sample3

Any observed barcode that does not match a barcode listed in the file (with a possible mismatch of 1) will be written to the sample specified by 'default'. Note that all barcodes must be of the same length.
""")
    parser.add_argument('output', metavar='output_basename', help="Base output directory, which must exist. As an example, if a given PE read has the barcode ACTACT, then it will be written to /output_basename/Sample_R1.fastq.gz and /output_basename/Sample_R2.fastq.gz")
    parser.add_argument('input', metavar='sample.fastq.gz', nargs="+", help="One or two gzipped fastq files. If two files are specified, it's assumed that the dataset is paired-end.")
    args = parser.parse_args(args)

    if args.sampleTable is None or args.output is None or args.input is None:
        parser.print_help()

    return args


def readSampleTable(sampleTable):
    """
    Read a sample table into a dict (keys are barcodes).
    Return the resulting dict and the barcode length.
    """
    f = open(sampleTable)
    d = dict()
    bcLen = 0
    for line in f:
        _ = line.strip().split("\t")
        d[_[0]] = _[1]
        if _[0] != 'default' and len(_[0]) > bcLen:
            bcLen = len(_[0])
    f.close()

    return (d, bcLen)


def matchSample(sequence, oDict, bcLen, buffer):
    """
    Match a barcode against the sample sheet with a possible edit distance of 1.

    Returns a tuple of the list of file pointers and True/False (whether the "default" output is being used)
    """
    bc = sequence[buffer:bcLen]
    if bc in oDict:
        return (bc, True)

    # Look for a 1 base mismatch
    for k, v in oDict.items():
        if ed.eval(k, bc) == 1:
            return (k, True)

    return ("default", False)


def writeRead(lineList, of, bc, args, doTrim=True):
    """
    Trim the read as specified and write to the appropriate file.

    Return the modified read name (for read #2, if needed)
    """
    rname = lineList[0]
    bcLen = len(bc)
    if doTrim:
        # Trim off the barcode
        lineList[1] = lineList[1][bcLen + args.buffer:]
        lineList[3] = lineList[3][bcLen + args.buffer:]

        # Fix the read name
        rname = rname.split()
        rname[0] = "{}_{}".format(rname[0], bc)
        rname = " ".join(rname)

    if rname[-1] != '\n':
        rname += '\n'

    of.write(rname)
    of.write(lineList[1])
    of.write(lineList[2])
    of.write(lineList[3])

    return rname


def writeRead2(lineList, of):
    # Fix the read name so it's read #2 rather than #1
    rname = lineList[0]
    rname = rname.split()
    rname[1] = "2{}".format(rname[1][1:])
    rname = " ".join(rname)

    if rname[-1] != '\n':
        rname += '\n'

    of.write(rname)
    of.write(lineList[1])
    of.write(lineList[2])
    of.write(lineList[3])


def processPaired(args, sDict, bcLen):
    f1 = gzip.GzipFile(args.input[0])
    f2 = gzip.GzipFile(args.input[1])

    for line1_1 in f1:
        line1_2 = f1.readline()
        line1_3 = f1.readline()
        line1_4 = f1.readline()
        line2_1 = f2.readline()
        line2_2 = f2.readline()
        line2_3 = f2.readline()
        line2_4 = f2.readline()

        (bc, isDefault) = matchSample(line1_2, sDict, bcLen, args.buffer)

        rname = writeRead([line1_1, line1_2, line1_3, line1_4], sDict[bc][0], bc, args, isDefault)
        writeRead2([rname , line2_2, line2_3, line2_4], sDict[bc][1])

    f1.close()
    f2.close()


def processSingle(args, sDict, bcLen):
    f1 = gzip.GzipFile(args.input[0])

    for line1_1 in f1.readline():
        line1_2 = f1.readline()
        line1_3 = f1.readline()
        line1_4 = f1.readline()
        (bc, isDefault) = matchSample(line1_2, sDict, bcLen, args.buffer)
        writeRead([line1_1, line1_2, line1_3, line1_4], sDict[bc], bc, args, isDefault)

    f1.close()


def main(args=None):
    args = parseArgs(args)

    sDict, bcLen = readSampleTable(args.sampleTable)

    # Open the output files
    if len(args.input) == 2:
        # of = [gzip.GzipFile("{}_R1.fastq.gz".format(args.output), "w"),
        #       gzip.GzipFile("{}_R2.fastq.gz".format(args.output), "w")]
        for k, v in sDict.items():
            sDict[k] = [gzip.GzipFile("{}/{}_R1.fastq.gz".format(args.output, v), "w"),
                        gzip.GzipFile("{}/{}_R2.fastq.gz".format(args.output, v), "w")]
    else:
        # of = [gzip.GzipFile("{}_R1.fastq.gz".format(args.output), "w")]
        for k, v in sDict.items():
            sDict[k] = [gzip.GzipFile("{}/{}_R1.fastq.gz".format(args.output, v), "w")]

    if len(args.input) == 2:
        processPaired(args, sDict, bcLen)
    else:
        processSingle(args, sDict, bcLen)


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["-h"]
    main(args)
