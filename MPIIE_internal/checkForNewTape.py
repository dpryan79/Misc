#!/usr/bin/env python
import argparse
import os
from math import ceil
from os.path import getsize, join
import glob
import sys

def getDSize(d) :
    """
    Returns the number of megabytes occupied by a directory tree. 

    Note that the number of blocks that a file occupies on the disk is used.
    """
    blockSize = os.statvfs("/data/rapidus").f_bsize
    o = 0
    for root, dirs, files in os.walk(d) :
        o += sum(ceil(getsize(join(root,name))/(blockSize)) for name in files)
        o += sum(ceil(getsize(join(root,dname))/(blockSize)) for dname in dirs)
    o += ceil(getsize(root)/(blockSize))
    o = o*blockSize/(1024*1024)
    return int(o)

def writeFiles(d) :
    """
    Write the *filelist and *dirlisting files.

    The *dirlisting file should be stored to tape since it contains the file sizes
    and permission of everything that will get backed up. The *filelist file
    contains a list of directories to back up (and the *dirlisting file).
    """
    if(len(d) > 1) :
        f = d[0].split("/")[-1].split("_")[2]
        t = d[-1].split("/")[-1].split("_")[2]
        odirs = "backup_to_tape_%s_through_%s.filelist" % (f, t)
        ofiles = "backup_to_tape_%s_through_%s.dirlisting" % (f, t)
        v = "/tmp/backup_to_tape_%s_through_%s.verbose" % (f, t)
        serr = "/tmp/backup_to_tape_%s_through_%s.stderr" % (f, t)
    elif(len(d) == 1) :
        f = d[0].split("/")[-1].split("_")[2]
        odirs = "backup_to_tape_%s.filelist" % f
        ofiles = "backup_to_tape_%s.dirlist" % f
        v = "/tmp/backup_to_tape_%s.verbose" % f
        serr = "backup_to_tape_%s.stderr" % f
    else :
        sys.stderr.write("Got an empty file list...which should only happen if a single directory is larger than 750 gigs!\n")
        return

    #This is fed to tar
    o = open(odirs, "w")
    o.write("%s\n" % ofiles)
    for dname in d :
        o.write("%s\n" % dname.split("/")[-1])
    o.close()

    #This is a list of files and their properties
    o = open(ofiles, "w")
    o.write("============= file contents of the following dirs on solserv2:/dont_touch_this/solexa_data/backup ====================\n")
    for dname in d :
        o.write("%s\n" % dname.split("/")[-1])
    o.write("================================================================================================================\n")
    o.close()
    for dname in d :
        cmd = "find %s -ls >> %s" % (dname, ofiles)
        os.system(cmd)

    #Print the command to be run to actually back things up
    sys.stdout.write("run:# mt -f /dev/st0 rewind; tar -cvp  -f /dev/nst0 --index-file=%s --files=%s 2> %s\n" % (v, odirs, serr))

def mainFunc(lastHiSeq, lastNextSeq) :
    """
    The main driver function.

    Receives the last directory written to tape, so it'll only check subsequent runs
    """
    maxSize = 798720 #780 gigs in megs
    curSize = 0

    # HiSeq
    dirs = glob.glob("/data/rapidus/sequencing_data/*_SN7001180_*")
    dirs.sort()
    dirsToBackup = []
    for d in dirs :
        runNum = int(d.split("/")[-1].split("_")[2])
        if(runNum <= lastHiSeq) :
            continue
        sz = getDSize(d)
        if(curSize + sz >= maxSize) :
            writeFiles(dirsToBackup)
            sys.stdout.write("This tape would contain %i megs of data\n" % curSize)
            curSize = sz
            dirsToBackup = [d]
        else :
            curSize += sz
            dirsToBackup.append(d)
    # NextSeq
    dirs = glob.glob("/data/rapidus/sequencing_data/*_NB501361_*")
    dirs.sort()
    for d in dirs :
        runNum = int(d.split("/")[-1].split("_")[2])
        if(runNum <= lastNextSeq ) :
            continue
        sz = getDSize(d)
        if(curSize + sz >= maxSize) :
            writeFiles(dirsToBackup)
            sys.stdout.write("This tape would contain %i megs of data\n" % curSize)
            curSize = sz
            dirsToBackup = [d]
        else :
            curSize += sz
            dirsToBackup.append(d)
    print("The remaining files occupy %i megs of data." % curSize)

parser = argparse.ArgumentParser(description="See if we have enough new data to write a tape.")
parser.add_argument('lastHiSeq', metavar='INT', type=int, help="The run number, for example 210 or 199, of the last backed up HiSeq flow cell.")
parser.add_argument('lastNextSeq', metavar='INT', type=int, help="The run number, for example 1 or 10, of the last backed up NextSeq flow cell.")
args = parser.parse_args()

if(args.lastHiSeq is None or args.lastNextSeq is None) :
    parser.print_help()
    sys.exit()

mainFunc(args.lastHiSeq, args.lastNextSeq)
