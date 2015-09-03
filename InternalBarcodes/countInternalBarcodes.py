#!/usr/bin/env python
import gzip
import sys

#CAAGCAGAAGACGGCATACGAGAT
#Added GTCCGC
indices = ["ATCACG","CGATGT", "TTAGGC", "TGACCA", "ACAGTG", "GCCAAT", "CAGATC", "ACTTGA", "GATCAG", "TAGCTT", "GGCTAC", "CTTGTA", "AGTCAA", "AGTTCC", "ATGTCA", "CCGTCC","GTCCGC", "GTGAAA", "GTGGCC", "GTTTCG", "CGTACG", "GAGTGG", "ACTGAT", "ATTCCT"]

def RC(idxs) :
    o = []
    for idx in idxs :
        idx = idx[::-1]
        idx = idx.replace("A","Q")
        idx = idx.replace("T","A")
        idx = idx.replace("Q","T")
        idx = idx.replace("C","Q")
        idx = idx.replace("G","C")
        idx = idx.replace("Q","G")
        o.append(idx)
    return o

indicesRev = RC(indices)
vals = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

def findBC(seq) :
    for i in range(len(indices)) :
        idx = seq.find("GAACTCCAGTCAC%s" % (indices[i]))
        if(idx >= 0) :
            vals[i] += 1
            return
    for i in range(len(indicesRev)) :
        idx = seq.find("%sGTGACTGGAGTTC" % (indicesRev[i]))
        if(idx >= 0) :
            vals[i] += 1
            return

f = gzip.open(sys.argv[1], "rb")
while 1:
    try :
        line1 = f.next()
        line2 = f.next()
        line3 = f.next()
        line4 = f.next()
    except :
        break

    findBC(line2)

for i in range(len(indices)) :
    print("%s\t%i" % (indices[i], vals[i]))
