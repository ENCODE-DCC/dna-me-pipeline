#!/bin/bash -e

if [ $# -ne 4 ]; then
    echo "usage v1: dname_bed_coor.sh <A_bed.gz> <B_bed.gz> <methyl_type> <out_root>"
    echo "Correlates two bed files containing methylation results.  Is independent of DX and ENCODE."
    echo "Requires bedtools intersectBed and path."
    # https://github.com/arq5x/bedtools2/releases
    exit -1; 
fi
A_bed_gz=$1          # First compressed bed file for correlation
B_bed_gz=$2          # Second compressed bed file for correlation
methyl_type=$3       # Type of methylation: CpG, CHG, CHH
out_root="${4}_${3}_corr" # The root of the output file (e.g. if "out" and type is "CpG" then "out_CpG_corr.txt" will be written).

echo "-- Expect to create '${out_root}.txt'"

echo "-- Correlate two bedmethyl files..."
set -x
intersectBed -a $A_bed_gz -b $B_bed_gz -wo | bedmethyl_corr.py - $methyl_type > ${out_root}.txt
set +x

echo "-- The results..."
#ls -l ${out_root}.txt
cat ${out_root}.txt

