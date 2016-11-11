#!/bin/bash -e

if [ $# -ne 1 ]; then
    echo "usage v1: bam-stats.sh <bam_file"
    echo "Generate 'samtools flagstat' and 'samtools stats' on a bam file."
    echo "Requires samtools on path."
    exit -1; 
fi
bam_file=$1  # Bam file to generate 'samtools flagstat' and 'samtools stats' on
bam_root=${bam_file%.bam}

echo "-- Expect to create '${bam_root}_flagstat.txt', '${bam_root}_stats.txt' and '${bam_root}_samstats_summary.txt'"

echo "-- Collect bam stats..."
set -x
samtools flagstat ${bam_root}.bam > ${bam_root}_flagstat.txt
samtools stats ${bam_root}.bam > ${bam_root}_samstats.txt
head -3 ${bam_root}_samstats.txt
grep ^SN ${bam_root}_samstats.txt | cut -f 2- > ${bam_root}_samstats_summary.txt
set +x
fi

echo "-- The results..."
ls -l ${bam_root}*
df -k .
