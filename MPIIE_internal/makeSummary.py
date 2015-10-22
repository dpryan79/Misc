#!/usr/bin/env python
import argparse
import ConfigParser
import os
import time
import csv
import subprocess
import glob
import xml.etree.ElementTree as et
import zipfile

def parseSampleSheet(fname) :
    o = dict()
    ename = ""
    e = []
    for line in open("/dont_touch_this/solexa_runs/%s/SampleSheet.csv" % fname.split("/")[-1]) :
        if(line.startswith(",")) :
            continue
        line = line.strip()
        if(line.startswith("[")) :
            if(len(e)) :
                o[ename] = e
            ename = line.split(",")[0]
            e = []
        else :
            e.append(line.split(","))
    if(len(e)>0):
        o[ename] = e
    if(len(o)):
        return o
    return None

def parseXML(dName) :
    try :
        t = et.parse("%s/Stats/DemultiplexingStats.xml" % dName).getroot()[0] #The first flow cell
    except :
        return None
    rv = []
    for proj in t.findall("Project") :
        if(proj.get("name") in ["default", "all"]) :
            continue
        for sam in proj.findall("Sample") :
            pname = proj.get("name")
            if(sam.get("name") in ["all"]) :
                continue
            sname = sam.get("name")
            for bc in sam.findall("Barcode") :
                if(bc.get("name") == "all") :
                    continue
                bname = bc.get("name")
                val = 0
                for lane in bc.findall("Lane") :
                    val += int(lane[0].text)
                rv.append([pname, sname, bname, val])
    if(len(rv) > 0) :
        return rv
    return None

def parseZipFiles(dname) :
    files = glob.glob("%s/FASTQC*/Sample*/*R1_fastqc.zip" % dname)
    rv = [] #[Project, library, %Duplication]
    for f in files :
        libName = f.split("/")[-2][7:]
        proj = f.split("/")[-3][15:]
        try :
            o = zipfile.ZipFile(f, "r")
        except :
            continue
        l = o.namelist()
        d = None
        for n in l :
            if(n.endswith("fastqc_data.txt")) :
                d = o.open(n)
                break
        if n is None :
            continue
        for line in d :
            line = line.strip().split("\t")
            if(line[0] == "#Total Deduplicated Percentage") :
                rv.append([proj, libName, 100-float(line[1])])
        d.close()
        o.close()
    if(len(rv)) :
        return rv
    return None

def parseFastqScreen(dname) :
    files = glob.glob("%s/Project*/Sample*/*R1_screen.txt" % dname)
    rv = [] #[Project, library, %human, %mouse, %fly, %fish, %adapter, %off-target (presumed)]
    for f in files :
        libName = f.split("/")[-2][7:]
        proj = f.split("/")[-3][8:]
        vals = [proj, libName, "0.0", "0.0", "0.0", "0.0", "0.0", "0.0"]
        pers = []
        for line in open(f) :
            if(line.startswith("#") or line.startswith("Library") or len(line) < 40) :
                continue
            line = line.strip().split("\t")
            if(line[0] == "Human") :
                vals[2] = line[5]
            elif(line[0] == "Mouse") :
                vals[3] = line[5]
            elif(line[0] == "Drosophila") :
                vals[4] = line[5]
            elif(line[0] == "Zebrafish") :
                vals[5] = line[5]
            elif(line[0] == "Adapter") :
                vals[6] = line[5]
            pers.append(float(line[5]))
        m = pers.index(max(pers))
        del pers[m]
        tot = 0.0
        for v in pers :
            tot += v
        vals[7] = "%f" % tot
        rv.append(vals)
    if(len(rv)) :
        return rv
    return None

def getIndices(header) :
    idxs = [None, None, None, None, None, None, None, None, None]
    targets = ["'ID'", "'User'", "'Group'", "'Organization'", "'Sample_Type'", "'Instrument'", "'Date_ProjectCreated'", "'Date_Submission'", "'Date_Finished'"]
    for i,v in enumerate(header) :
        if(v in targets) :
            idxs[targets.index(v)] = i
    return idxs

def subsetLine(line, indices) :
    o = []
    for i in indices :
        o.append(line[i].replace("'",""))
    if(len(o)) :
        return o
    return None

#Projects can have multiple flow cells, but apparently not multiple sets of dates?!?!?
def parseMDBtools() :
    rv = []
    fp = subprocess.Popen(["/home/pipegrp/bin/mdb-export","-Q", "-d","'\t'","/eva_data/SAMBA/seq_share/01_MasterTool/ProjectDB.accdb","Data"], stdout = subprocess.PIPE).stdout
    indices = getIndices(fp.next().split("\t"))
    for line in fp :
        line = line.split("\t")
        if(line[0] == "ID_no") :
            indices = getIndices(line)
        rv.append(subsetLine(line, indices))
        if(rv[-1] is None) :
            del rv[-1]
    if(len(rv)) :
        return rv
    return None

class Sample:
    '''
    I hold a single sample's information
    '''
    def __init__(self) :
        self.LibraryID = ""
        self.SampleName = ""
        self.Index = ""
        self.nFragments = "0"
        self.PercentDuplication = "0.0"
        self.PercentSpecies = ["0.0", "0.0", "0.0", "0.0", "0.0", "0.0"] #human, mouse, fly, zebrafish, adapter, off_species

    def title(self) :
        rv = ["LibraryID", "SampleName"]
        rv.extend(["Index", "#Fragments","PercentDuplication","PercentHuman", "PercentMouse", "PercentFly", "PercentZebraFish", "PercentAdapter", "PercentOffSpecies"])
        return rv

    def __repr__(self) :
        rv = [self.LibraryID, self.SampleName]
        rv.extend([self.Index, self.nFragments, self.PercentDuplication])
        rv.extend(self.PercentSpecies)
        return rv

class Project:
    '''
    I'm a project, I hold multiple samples
    '''
    def __init__(self):
        self.ProjectID = ""
        self.Submitter = ""
        self.PI = ""
        self.Organization = ""
        self.DateSubmitted = "unset" #should be per sample!
        self.DatePrepped = "unset" #should be per sample! #mdbtools?
        self.TargetSpecies = "unset" #mdbtools?
        self.LibraryType = "" #should be per sample!
        self.FragmentLength = "unset" #should be per sample! #mdbtools?
        self.RunType = ""
        self.Samples = []

    def title(self) :
        rv = (["ProjectID", "Submitter", "PI", "Organization"])
        rv.extend(["DateSubmitted", "DatePrepared", "TargetSpecies", "LibraryType", "FragmentLength", "RunType"])
        rv.extend(Sample().title())
        return rv

    def __repr__(self) :
        o = []
        if(len(self.Samples)) :
            for i in range(len(self.Samples)) :
                rv = [self.ProjectID, self.Submitter, self.PI, self.Organization]
                rv.extend([self.DateSubmitted, self.DatePrepped, self.TargetSpecies, self.LibraryType, self.FragmentLength, self.RunType])
                rv.extend(self.Samples[i].__repr__())
                o.append(rv)
            return o
        else :
            return [self.ProjectID, self.Submitter, self.PI, self.Organization]

    def sampleIndex(self, sName, idx) :
        for i in range(len(self.Samples)) :
            if(sName == self.Samples[i].SampleName) :
                if(idx == self.Samples[i].Index) :
                    return i
        return None

    def libIndex(self, libName, idx) :
        for i in range(len(self.Samples)) :
            if(libName == self.Samples[i].LibraryID) :
                if(idx == self.Samples[i].Index) :
                    return i
        return None

    #This is just libIndex without the barcode
    def libIndex2(self, libName) :
        for i in range(len(self.Samples)) :
            if(libName == self.Samples[i].LibraryID) :
                return i

    def addSample(self, line) :
        s = Sample()
        s.LibraryID = line[1]
        s.SampleName = line[2]
        s.Index = line[6]
        self.Samples.append(s)
        return len(self.Samples)-1

class FlowCell:
    '''
    I'm a flow cell, I hold multiple projects
    '''
    def __init__(self):
        self.RunID = ""
        self.FlowCellID = ""
        self.Instrument = ""
        self.DateSequenced = ""
        self.DateDelivered = ""
        self.ReadLength1 = "0"
        self.ReadLength2 = "0"
        self.Projects = []

    def title(self) :
        rv = ["RunID", "FlowCellID", "Instrument", "DateSequenced", "DateDelivered", "ReadLength1", "ReadLength2"]
        rv.extend(Project().title())
        return "\t".join(rv)
        

    def __repr__(self) :
        if(len(self.Projects) > 0) :
            o = []
            for i in range(len(self.Projects)) :
                tempO = self.Projects[i].__repr__()
                for e in tempO :
                    foo = [self.RunID, self.FlowCellID, self.Instrument, self.DateSequenced, self.DateDelivered, self.ReadLength1, self.ReadLength2]
                    foo.extend(e)
                    o.append("\t".join(foo))
        else :
            o = [self.RunID, self.FlowCellID, self.Instrument, self.DateSequenced, self.DateDelivered, self.ReadLength1, self.ReadLength2, self.Projects]
            o = "\t".join(o)
        return "\n".join(o)

    def getStats(self, fname):
        self.RunID= fname.split("_")[4]
        self.FlowCellID = fname.split("_")[-1]
        self.Instrument = fname.split("_")[3]
        self.DateSequenced = fname.split("_")[2].split("/")[1]
        self.DateDelivered = time.strftime("%y%m%d", time.gmtime(os.path.getctime(fname))) #Creation date

    def projectIndex(self, projName) :
        for i in range(len(self.Projects)) :
            if(projName == self.Projects[i].ProjectID) :
                return i
        return None

    def hasSample(self, projName, libName, Index) :
        idx = self.projectIndex(projName)
        if(idx is not None) :
            if(self.Projects[idx].libIndex(libName, Index) is not None) :
                return True
            return False
        return False

    def addProject(self, line) :
        p = Project()
        if(len(line) == 9) :
            p.ProjectID = line[7]
        elif(len(line) == 7) :
            p.ProjectID = line[5]
        else :
            assert(1==0)
        p.addSample(line)
        self.Projects.append(p)
        return len(self.Projects)-1

    def addSample(self, line) :
        if(len(line) == 9) :
            idx = self.projectIndex(line[7])
        elif(len(line) == 7) :
            idx = self.projectIndex(line[5])
        else :
            assert(1==0)

        if(idx is None) :
            idx = self.addProject(line)
        else :
            self.Projects[idx].addSample(line)

    #TODO: Handle no sample sheet (i.e., some rapid runs)
    def addSampleSheet(self, ss) :
        if "[Reads]" in ss :
            reads = ss["[Reads]"]
            self.ReadLength1 = reads[0][0]
            if(len(reads) > 1) :
                self.ReadLength2 = reads[1][0]
        if "[Data]" in ss :
            d = ss["[Data]"]
            for line in d:
                #Do we already have this sample?
                if(line[0] == "Lane"):
                    continue
                if(len(line) == 9) :
                    if(self.hasSample(line[7], line[1], line[6]) == False) :
                        self.addSample(line)
                elif(len(line) == 7) :
                    if(self.hasSample(line[5], line[1], "NA") == False) :
                        self.addSample(line)

    def setPercentDupes(self, x) :
        for line in x :
            i = self.projectIndex(line[0])
            if i is None :
                continue
            j = self.Projects[i].libIndex2(line[1])
            if j is None :
                continue
            self.Projects[i].Samples[j].PercentDuplication = "%f" % (line[2])

    def setContamination(self, x) :
        for line in x :
            i = self.projectIndex(line[0])
            if i is None :
                continue
            j = self.Projects[i].libIndex2(line[1])
            if j is None :
                continue
            self.Projects[i].Samples[j].PercentSpecies = line[2:]

    def setMDBvalues(self, x) :
        for line in x :
            i = self.projectIndex(line[0])
            if i is None :
                continue
            self.Projects[i].Submitter = line[1]
            self.Projects[i].PI = line[2]
            self.Projects[i].Organization = line[3]
            self.Projects[i].DateSubmitted = line[7]
            self.Projects[i].LibraryType = line[4]
            self.Projects[i].RunType = line[5]

    def setNFragments(self, x) :
        for e in x :
            i = self.projectIndex(e[0])
            if i is None :
                continue
            j = self.Projects[i].sampleIndex(e[1],e[2])
            if j is None :
                continue
            self.Projects[i].Samples[j].nFragments = "%i" % (e[3])

def getRuns(lastRun) :
    runs = glob.glob("/eva_data/sequencing_data/*SN7001180*")
    newRuns = []
    for run in runs :
        id = int(run.split("_")[4])
        if(id >= lastRun) :
            newRuns.append(run)
    return newRuns

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--first", type=int, help="The run number from which you want to produce output (e.g., 232). The default is 0 (include all runs).", default=0)
args = parser.parse_args()

runs = getRuns(args.first)
if(len(runs) > 0) :
    print(FlowCell().title())
runs.sort()

#Parse mdbtools' dump of the Access database
MDBvalues = parseMDBtools()

for run in runs :
    FC = FlowCell()

    #Get the file name and fstat stats
    FC.getStats(run)

    #Parse the sample sheet (producing the projects and samples
    ss = parseSampleSheet(run)
    if(ss) :
        FC.addSampleSheet(ss)

    #Parse the XML stats
    x = parseXML(run)
    if(x) :
        FC.setNFragments(x)

    #Parse FastQC
    stats = parseZipFiles(run)
    if(stats) :
        FC.setPercentDupes(stats)

    #Parse fastq_screen
    stats = parseFastqScreen(run)
    if(stats) :
        FC.setContamination(stats)

    #Add the database info
    if(MDBvalues) :
        FC.setMDBvalues(MDBvalues)

    #Append the output to a file
    print(FC)
