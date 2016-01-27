#!/usr/bin/env python
import argparse
import subprocess
import glob
import xml.etree.ElementTree as et

def parseXML(dName) :
    try :
        t = et.parse("%s/Stats/DemultiplexingStats.xml" % dName).getroot()[0] #The first flow cell
    except :
        return None

    rv = []
    for proj in t.findall("Project") :
        if(proj.get("name") == "all") :
            # sample
            for sam in proj.findall("Sample"):
                if(sam.get("name") != "all"):
                    continue
                # barcode
                for bc in sam.findall("Barcode") :
                    if(bc.get("name") != "all") :
                        continue
                    for lane in bc.findall("Lane") :
                        rv.append([lane.get("number"), lane[0].text])
    if(len(rv) > 0) :
        return rv
    return None

def getRuns(lastRunHiSeq, lastRunNextSeq) :
    runs = glob.glob("/eva_data/sequencing_data/*SN7001180*")
    newRuns = []
    for run in runs :
        id = int(run.split("_")[4])
        if(id >= lastRunHiSeq) :
            newRuns.append(run)

    runs = glob.glob("/eva_data/sequencing_data/*NB501361*")
    for run in runs :
        id = int(run.split("_")[4])
        if(id >= lastRunNextSeq) :
            newRuns.append(run)

    return newRuns


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--firstHiSeq", type=int, help="The HiSeq run number from which you want to produce output (e.g., 232). The default is 0 (include all runs).", default=0)
parser.add_argument("-F", "--firstNextSeq", type=int, help="The NextSeq run number from which you want to produce output (e.g., 232). The default is 0 (include all runs).", default=0)
args = parser.parse_args()

runs = getRuns(args.firstHiSeq, args.firstNextSeq)
runs.sort()

stats = []
for run in runs :
    o = parseXML(run)
    if o is not None:
        stats.append([run.split("/")[3], o])

if len(stats) > 0:
    print("Flow cell\tMachine\tLane\tPassing cluster count")
    for l in stats:
        machine = "HiSeq"
        if l[0][7:9] == "NB":
            machine = "NextSeq"
        if len(l[1]) == 2:
            machine += " Rapid Run"
        FC = l[0]
        for lane in l[1]:
            print("{}\t{}\t{}\t{}".format(FC, machine, lane[0], lane[1]))
