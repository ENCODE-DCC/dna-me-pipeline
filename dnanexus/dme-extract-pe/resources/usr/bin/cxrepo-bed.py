#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2014 Junko Tsuji

# This script generates bed files as final outputs of
# DNA methylation pipeline by taking a CX_report file
# generated from bismark_methylation_extractor.

from optparse import OptionParser
import sys, os.path, string

# RGB strings and integer values
rgbDict = { (0 , 5) :   "0,255,0", # 65280,
            (6 , 15):  "55,255,0", # 3669760,
            (16, 25): "105,255,0", # 6946560,
            (26, 35): "155,255,0", # 10223360,
            (36, 45): "205,255,0", # 13500160,
            (46, 55): "255,255,0", # 16776960
            (56, 65): "255,205,0", # 16764160,
            (66, 75): "255,155,0", # 16751360,
            (76, 85): "255,105,0", # 16738560,
            (86, 95):  "255,55,0", # 16725760,
            (96,100):   "255,0,0"  # 16711680
           }

def write(w, line, name):
    chrom      = line[0]
    chromStart = line[1]
    chromEnd   = str( int(line[1])+1 )
    strand = line[2]
    meth   = float(line[3])
    nometh = int(line[4])
    readCount = int(meth) + nometh
    itemRgb = rgbDict[(0,5)]
    if (meth + nometh) > 0:
        percentMeth = int(round(meth/(meth+nometh)*100))
        for L in rgbDict:
            if L[0] <= percentMeth and percentMeth <= L[1]:
                itemRgb = rgbDict[L]
                break
    else:
        percentMeth = 0
    w.write("\t".join([ chrom,
                        chromStart,
                        chromEnd,
                        name,
                        str(min(readCount, 1000)),
                        strand,
                        chromStart,
                        chromEnd,
                        itemRgb,
                        str(readCount),
                        str(percentMeth) ])+"\n")


def outAll(files, f, name):
    for line in open(f):
        line = line.rstrip("\n").split("\t")
        write(files[line[5]], line, name)


def outExact(files, f, name):
    for line in open(f):
        line = line.rstrip("\n").split("\t")
        if   line[6][1] == "N":
            continue
        elif line[5] != "CG":
            if line[6][2] == "N":
                continue
        write(files[line[5]], line, name)


def cxrepoBed(opts, args):
    f = os.path.basename(args[0])
    files = {'CG' :open(opts.output+"/CG_" +f, "w"),
             'CHG':open(opts.output+"/CHG_"+f, "w"),
             'CHH':open(opts.output+"/CHH_"+f, "w")}
    if opts.ns == True:
        outAll(files, args[0], opts.name)
    else:
        outExact(files, args[0], opts.name)
    for of in files:
        files[of].close()


if __name__ == "__main__":
    prog = os.path.basename(sys.argv[0])
    fPrefix = os.path.basename(sys.argv[-1]).split(".CX_report")[0]
    usage = "%prog [option] cxReport"
    description = "Convert a 'CX_report' file to bed files"

    op = OptionParser(usage=usage, description=description)

    op.add_option("-n","--output-Ns",
                  dest="ns",
                  action="store_true",
                  default=False,
                  help="Output cytosine contexts in which include 'N' (default=%default)")

    op.add_option("-N","--bedmethyl-name",
                  dest="name",
                  type="string",
                  action="store",
                  default=fPrefix,
                  help="Insert name into the 3rd column of the bedMethyl format (default=[cxReport filename])",
                  metavar="NAME")

    op.add_option("-o","--output-place",
                  dest="output",
                  type="string",
                  action="store",
                  default="./",
                  help="Place for output files (default='./')",
                  metavar="PATH")

    (opts, args) = op.parse_args()

    try:
        cxrepoBed(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        sys.exit(prog + ": error: " + str(e))
