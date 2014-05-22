#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2014 Junko Tsuji

# This script generates bed files as final outputs of
# DNA methylation pipeline by taking a CX_report file
# generated from bismark_methylation_extractor.

from optparse import OptionParser
import sys, os.path, string

def write(w, line):
    chrom  = line[0]
    beg    = line[1]
    end    = str( int(line[1])+1 )
    strand = line[2]
    sum    = float(line[3])+float(line[4])
    genome = line[6]
    if sum > 0:
        mc = "{0:.4f}".format(float(line[3])/sum*100)
    else:
        mc = "0.0000"
    w.write("\t".join([ chrom, beg, end, genome, mc, strand ])+"\n")


def outAll(files, f):
    for line in open(f):
        line = line.rstrip("\n").split("\t")
        write(files[line[5]], line)


def outExact(files, f):
    for line in open(f):
        line = line.rstrip("\n").split("\t")
        if   line[6][1] == "N":
            continue
        elif line[5] != "CG":
            if line[6][2] == "N":
                continue
        write(files[line[5]], line)


def cxrepoBed(opts, args):
    f = os.path.basename(args[0])
    files = {'CG' :open(opts.output+"/CG_" +f, "w"),
             'CHG':open(opts.output+"/CHG_"+f, "w"),
             'CHH':open(opts.output+"/CHH_"+f, "w")}
    if opts.ns == True:
        outAll(files, args[0])
    else:
        outExact(files, args[0])
    for of in files:
        files[of].close()


if __name__ == "__main__":
    prog = os.path.basename(sys.argv[0])
    usage = "%prog [option] cxReport"
    description = "Convert a 'CX_report' file to bed files"

    op = OptionParser(usage=usage, description=description)

    op.add_option("-n","--output-Ns",
                  dest="ns",
                  action="store_true",
                  default=False,
                  help="Output cytosine contexts in which include 'N' (default=%default)")

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
