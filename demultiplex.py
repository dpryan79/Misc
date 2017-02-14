#!/usr/bin/env python
import gzip
import os
import editdistance

files = dict()
barcodes = dict()

def fillLines(f1, I1, I2, f2):
    out1 = []
    out2 = []
    outI1 = None
    outI2 = None
    for _ in range(4):
        out1.append(f1.next())
        out2.append(f2.next())
        if _ == 1:
            outI1 = I1.next()
            outI2 = I2.next()
        else:
            I1.next()
            I2.next()
    return (out1, outI1, outI2, out2)


def getSampleName(i1, i2, barcodes):
    # Try a 6 base barcode
    if i1[:6] in barcodes:
        return barcodes[i1[:6]]

    if i1 in barcodes:
        if i2 in barcodes[i1]:
            return barcodes[i1][i2]

    for k, v in barcodes.items():
        if len(k) == 6:
            if editdistance.eval(k, i1[:6]) == 1:
                return v
        else:
            if editdistance.eval(k, i1) <= 2:
                for k2, v2 in v.items():
                    if editdistance.eval(k2, i2) <= 2:
                        return v2
    return None


for line in open("samples.txt"):
    _ = line.strip().split(",")
    odir = "Project_{}/Sample_{}".format(_[9], _[1])
    try:
        os.makedirs(odir)
    except:
        pass
    of1 = gzip.GzipFile("{}/{}_R1.fastq.gz".format(odir, _[2]), "w")
    of2 = gzip.GzipFile("{}/{}_R2.fastq.gz".format(odir, _[2]), "w")

    files[_[2]] = [of1, of2]
    if len(_[6]) == 6:
        barcodes[_[6]] = _[2]
    else:
        if _[6] in barcodes:
            barcodes[_[6]][_[8]] = _[2]
        else:
            barcodes[_[6]] = dict()
            barcodes[_[6]][_[8]] = _[2]

in1 = gzip.GzipFile("Undetermined_S0_R1_001.fastq.gz")
in2 = gzip.GzipFile("Undetermined_S0_I1_001.fastq.gz")
in3 = gzip.GzipFile("Undetermined_S0_I2_001.fastq.gz")
in4 = gzip.GzipFile("Undetermined_S0_R2_001.fastq.gz")

unknown1 = gzip.GzipFile("Unknown_R1.fastq.gz", "w")
unknown2 = gzip.GzipFile("Unknown_R2.fastq.gz", "w")
while True:
    try:
        (lines1, i1, i2, lines2) = fillLines(in1, in2, in3, in4)
        sampleName = getSampleName(i1.strip(), i2.strip(), barcodes)
        if sampleName is not None:
            files[sampleName][0].write("{}{}{}{}".format(*lines1))
            files[sampleName][1].write("{}{}{}{}".format(*lines2))
        else:
            unknown1.write("{}{}{}{}".format(*lines1))
            unknown2.write("{}{}{}{}".format(*lines2))
    except:
        break

#Close things up
unknown1.close()
unknown2.close()
in1.close()
in2.close()
in3.close()
in4.close()
for k,v in files.items():
    v[0].close()
    v[1].close()
