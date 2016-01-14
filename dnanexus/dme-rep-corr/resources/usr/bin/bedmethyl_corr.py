#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script calculates pearson correlation coefficients
# of two bedMethyl files.

import sys, os.path, fileinput
from optparse import OptionParser
from math import sqrt

def bedmethylCorr(args):
    pSum = 0
    sum1, sum2 = 0, 0
    sqSum1, sqSum2 = 0, 0
    n, a = 0, 0
    for i, x in enumerate(fileinput.input(args[0])):
        x = x.rstrip().split("\t")
        a = i
        r1, r2 = int(x[4]), int(x[15])  # number of reads
        m1, m2 = float(x[10]), float(x[21]) # methylation level
        if r1 < 10 or r2 < 10: 
            continue
        sum1 += m1
        sum2 += m2
        sqSum1 += pow(m1,2)
        sqSum2 += pow(m2,2)
        pSum += m1*m2
        n += 1
    den = 0
    if n > 0:
        num = pSum - (sum1*sum2/n)
        den = sqrt((sqSum1 - pow(sum1,2)/n) * (sqSum2 - pow(sum2,2)/n))
    if den == 0: 
        return 0, n, a + 1
    return num/den, n, a + 1

if __name__ == "__main__":
    usage = "intersectBedStdin -a A.bed -b B.bed -wo | %prog - 'CpG'\nversion: v1"
    description = "Calculates pearson correlation coefficient"

    op = OptionParser(usage=usage, description=description)
    (opts, args) = op.parse_args()

    if len(args) != 2: 
        op.error("input should be piped with 'intersectBed -a bedMethyl1 -b bedMethyl -wo'")

    try:
        c, n, a = bedmethylCorr(args)
        print "Pearson Correlation Coefficient=%.4f\n%s pairs with atleast 10 reads each=%d\n%s pairs=%d" % \
                                                            (c, sys.argv[2], n, sys.argv[2], a)
    except KeyboardInterrupt: 
        pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))

