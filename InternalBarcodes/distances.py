#!/usr/bin/env python
import Levenshtein as L

indices = ["ATCACG","CGATGT", "TTAGGC", "TGACCA", "ACAGTG", "GCCAAT", "CAGATC", "ACTTGA", "GATCAG", "TAGCTT", "GGCTAC", "CTTGTA", "AGTCAA", "AGTTCC", "ATGTCA", "CCGTCC","GTCCGC", "GTGAAA", "GTGGCC", "GTTTCG", "CGTACG", "GAGTGG", "ACTGAT", "ATTCCT"]

for i in range(len(indices)-1) :
    for j in range(i+1,len(indices)) :
        d = L.distance(indices[i],indices[j])
        #if(d < 4) :
        #    print("%i\t%s\t%s" % (d,indices[i], indices[j]))
        print("%i\t%s\t%s" % (d,indices[i], indices[j]))
