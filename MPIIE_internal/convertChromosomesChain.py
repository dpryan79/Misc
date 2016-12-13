#!/usr/bin/env python
import csv
import argparse
import sys

def parseArguments(args=None) :
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Convert chromosome names between systems (e.g., UCSC to Ensembl) in a chain file. This program accepts a two column tab separated file as input. Column 1 of the file specifies the chromosome name of the input. Column 2 specifies what that name will be converted to. All entries in the file MUST be present or the program will crash. Lines starting with # are ignored. Output is to stdout. The file to be converted is accepted via a pipe.")
    parser.add_argument('mapFile', metavar='mapFile', help="The tab-separated file containing the chromosome name mappings in the first genome")
    parser.add_argument('mapFile2', metavar='mapFile2', help="The tab-separated file containing the chromosome name mappings in the second genome")
    args = parser.parse_args()
    if(args.mapFile == None) :
        parser.print_help()
        sys.exit()
    return(args)

def makeMappings(fname) :
    mappings = dict()
    try: 
        for line in csv.reader(open(fname, "r"), dialect="excel-tab") :
            assert line[0] not in mappings
            if line[1] != "" :
                mappings[line[0]] = line[1]
            else :
                continue
    except:
        sys.exit("Couldn't open %s!\n" % fname)
    return(mappings)

def main(args) :
    mappings = makeMappings(args.mapFile)
    mappings2 = makeMappings(args.mapFile2)
    out = sys.stdout
    lastKey = None
    lastVal = None
    for line in sys.stdin:
        if line.startswith("chain"):
            cols = line.strip().split()
            if cols[2] not in mappings:
                sys.exit("%s not in mappings!\n" % cols[2])
            cols[2] = mappings[cols[2]]
            if cols[7] not in mappings2:
                sys.exit("%s not in mappings!\n" % cols[7])
            cols[7] = mappings2[cols[7]]
            line = " ".join(cols)
            out.write("{}\n".format(line))
        else:
            out.write("{}".format(line))


if __name__ == "__main__" :
    args = parseArguments()
    main(args)
