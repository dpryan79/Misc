#!/usr/bin/env python
import gzip
import hashlib
import argparse
import sys

def getArgs() :
    parser = argparse.ArgumentParser(description="Compare 2 gzipped fastq files according to the md5 of their reads.")
    parser.add_argument("f1", metavar="fastq1.gz", help="The first gzipped fastq file.")
    parser.add_argument("f2", metavar="fastq2.gz", help="The second gzipped fastq file.")
    parser.add_argument("o1", metavar="different1.gz", help="Reads from fastq1.gz not present or different in fastq2.gz are written here.")
    parser.add_argument("o2", metavar="different2.gz", help="As above, but for reads from fastq2.gz.")
    parser.add_argument("-l", metavar="rlen", type=int, default=0, help="The maximum number of bases in each read to consider.")
    args = parser.parse_args()
    if(args.f1 is None or args.f2 is None or args.o1 is None or args.o2 is None) :
        args.print_help()
        sys.exit()
    return args

def getRead(fp, maxLen) :
    line = fp.next().split(" ")[0].split(":")
    line = line[1]+":"+line[2]+":"+line[3]+":"+line[4]+":"+line[5]+":"+line[6]
    line += "\t"
    line = ""
    if(maxLen > 0) :
        line += fp.next().strip()[:maxLen]
        line += "\t"
        fp.next()
        line += fp.next().strip()[:maxLen]
    else :
        line += fp.next().strip()
        line += "\t"
        fp.next()
        line += fp.next().strip()
    return line

def writeFQ(fp, read) :
    r = read.split("\t")
    fp.write("@read\n%s\n+\n%s\n" % (r[0],r[1]))

def main(args) :
    ht = set()
    
    #Read in file #1
    f1 = gzip.open(args.f1, "rb")
    while 1:
        try : 
            md5 = hashlib.md5()
            read = getRead(f1, args.l)
            md5.update(read)
            h = md5.digest()
            ht.add(h)
        except :
            break
    f1.close()
    print("Finished reading in fastq #1")

    #Compare file #2 to #1
    f2 = gzip.open(args.f2, "rb")
    o2 = gzip.open(args.o2, "wb")
    identical = 0
    different = 0
    while 1:
        try : 
            md5 = hashlib.md5()
            read = getRead(f2, args.l)
            md5.update(read)
            h = md5.digest()
            if h in ht :
                ht.remove(h)
            else :
                different += 1
                writeFQ(o2, read)
        except :
            break
    f2.close()
    o2.close()
    print("%s vs. %s\t Identical %i\tDifferent %i" % (args.f1, args.f2, identical, different))

    #Remaining reads in the hash table
    f1 = gzip.open(args.f1, "rb")
    o1 = gzip.open(args.o1, "wb")
    while 1:
        try :
            md5 = hashlib.md5()
            read = getRead(f1, args.l)
            md5.update(read)
            h = md5.digest()
            if h in ht :
                writeFQ(o1, read)
        except :
            break
    f1.close()
    o1.close()

if __name__ == "__main__" :
    args = getArgs()
    main(args)
