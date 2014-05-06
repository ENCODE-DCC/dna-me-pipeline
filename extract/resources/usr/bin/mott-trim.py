#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2014 Junko Tsuji

# This script computes scores based on running sum
# algorithm (mott trimming algorithm) to trim 3'
# ends of reads, with the score below threshould.

from optparse import OptionParser
import sys, os.path, re

nPat = re.compile('N', re.IGNORECASE)

def printit(entry, cutoff, base, minlen):
    index = -1
    minValue = 100000000
    scoreValue = 0
    rawseq = entry[1][::-1]
    rawqual = entry[3][::-1]
    for i in range(len(rawseq)):
        qscore = ord(rawqual[i]) - base
        scoreValue += (qscore - cutoff)
        if scoreValue < minValue:
            minValue = scoreValue
            index = i + 1
    if minValue > cutoff:
        seq = rawseq[::-1]
        qual = rawqual[::-1]
    else:
        seq = rawseq[index:][::-1]
        qual = rawqual[index:][::-1]
    L = len(seq)
    part = "0,"+str(L)
    if len(nPat.findall(seq)) != L and L >= minlen:
        print "\n".join([entry[0]+":"+part, seq, entry[2], qual])


def checkargs(opts):
    if opts.fqtype.lower() == "sanger":
        return 33
    if opts.fqtype.lower() == "illumina":
        return 64
    raise Exception("input correct fastq type")


def mottTrim(opts, args):
    entry = []
    lineCnt = 0
    base = checkargs(opts)
    for line in open(args[0]):
        line = line.rstrip("\n")
        if line == "":
            continue
        lineCnt += 1
        entry.append(line)
        if lineCnt == 4:
            printit(entry, opts.cutoff, base, opts.minlen)
            entry = []
            lineCnt = 0


if __name__ == "__main__":
    prog = os.path.basename(sys.argv[0])
    usage = "%prog [options] readFastq"
    description = "Trim 3' or both ends of reads using mott trimming algorithm"

    op = OptionParser(usage=usage, description=description)
    
    op.add_option("-q","--quality-cutoff",
                  dest="cutoff",
                  type="int",
                  action="store",
                  default=3,
                  help="PHRED quality score cutoff (default=%default)",
                  metavar="VALUE")
    op.add_option("-m", "--min-length",
                  dest="minlen",
                  type="int",
                  action="store",
                  default=30,
                  help="Minimum length of trimmed reads (default=%default)",
                  metavar="VALUE")
    op.add_option("-t", "--fastq-type",
                  dest="fqtype",
                  type="string",
                  action="store",
                  default="sanger",
                  help="PHRED quality score format: 'sanger', 'illumina' (default='%default')",
                  metavar="TYPE")

    (opts, args) = op.parse_args()

    try:
        mottTrim(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        sys.exit(prog + ": error: " + str(e))
