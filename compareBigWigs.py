#!/usr/bin/env python
import pyBigWig
import argparse

parser = argparse.ArgumentParser(description="Compare two bigWig files and ensure their intervals/headers are the same.")
parser.add_argument("bw", metavar="bigWig files", help="bigWig files to compare, there must be exactly 2 of them", nargs=2)
args = parser.parse_args()

bw1 = pyBigWig.open(args.bw[0])
bw2 = pyBigWig.open(args.bw[1])
assert(bw1 is not None)
assert(bw2 is not None)

#assert(bw1.header() == bw2.header()) #N.B., the zoomLevels need to be the same
#assert(bw1.chroms() == bw2.chroms())

for k,v in bw1.chroms().iteritems():
    ints1 = bw1.intervals(k)
    ints2 = bw2.intervals(k)
    assert(ints1 == ints2)

#To Do, check the zoom levels

bw1.close()
bw2.close()
