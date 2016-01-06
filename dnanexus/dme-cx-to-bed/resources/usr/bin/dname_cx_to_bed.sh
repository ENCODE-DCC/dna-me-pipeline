#!/bin/bash -e

if [ $# -ne 2 ]; then
    echo "usage v1: dname_cx_to_bed.sh <cx_report> <chrom.sizes>"
    echo "Converts coverage2cytosine CX_report to beds and bigBeds.  Is independent of DX and ENCODE."
    exit -1; 
fi
cx_report=$1   # uncompressed CX_report output of coverage2cytosine.
chrom_sizes=$2 # chrom.sizes file used by bedGraphToBigWig
target_root=${cx_report%.CX_report.txt}

echo "-- Create beds..."
set -x
mkdir -p output
mv $cx_report output/
cxrepo-bed.py -o output/ output/${cx_report}
mv output/CG_${target_root}.CX_report.txt  ${target_root}_CpG.bed
mv output/CHG_${target_root}.CX_report.txt ${target_root}_CHG.bed
mv output/CHH_${target_root}.CX_report.txt ${target_root}_CHH.bed
set +x

echo "-- Convert to bigBed..."
set -x
bedToBigBed ${target_root}_CHH.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 $chrom_sizes ${target_root}_CHH.bb
pigz ${target_root}_CHH.bed
bedToBigBed ${target_root}_CHG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 $chrom_sizes ${target_root}_CHG.bb
pigz ${target_root}_CHG.bed
bedToBigBed ${target_root}_CpG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 $chrom_sizes ${target_root}_CpG.bb
pigz ${target_root}_CpG.bed
set +x

echo "-- The results..."
ls -l ${target_root}*

