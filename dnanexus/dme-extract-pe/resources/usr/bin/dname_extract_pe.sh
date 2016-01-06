#!/bin/bash -e

if [ $# -lt 4 ] || [ $# -gt 4 ]; then
    echo "usage v1: dname_extract_pe.sh <index.tgz> <bismark.bam> <ncpus> [--scorched_earth]"
    echo "Extracts methylation from paired-end bismark bam.  Is independent of DX and ENCODE."
    echo "If --scorched_earth will remove everything, including input bam and index in order to maximize available storage."
    exit -1; 
fi
index_tgz=$1   # Index archive containing input/*.fa reference, input/Bisulfite_Genome/* index made with chosen bowtie version.
               # Also input/lambda/*, input/lambda/Bisulfite_Genome/*.  Will be expanded and REMOVED for storage efficiency
bismark_bam=$2 # bam input aligned with bismark and bowtie version that matches index_tgz.
ncpus=$3       # Number of cores available for bismark --multi
scorched="no"  # --scorched_earth will remove everything, including input bam and index in order to maximize available storage.
if [ $# -eq 4 ] && [ "$4" == "--scorched_earth" ]; then
    scorched="earth"
fi
target_root=${bismark_bam%.bam}

echo "-- Uncompress index archive..."
set -x
tar zxvf $index_tgz
rm $index_tgz
set +x
bowtie_ver="bowtie1"
if [ -f input/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2 ]; then
    bowtie_ver="bowtie2"
fi

echo "-- Collect bam stats..."
set -x
samtools flagstat ${target_root}.bam > ${target_root}_flagstat.txt
set +x
# NOTE: samtools stats may take longer than it is worth 
#samtools stats ${target_root}.bam > ${target_root}_samstats.txt
#head -3 ${target_root}_samstats.txt
#grep ^SN ${target_root}_samstats.txt | cut -f 2- > ${target_root}_samstats_summary.txt

# NOTE: Better to use sam and let extractor use more threads, but this takes up precious storage
alignments_file=${target_root}.bam
# What bam_size constitutes too large for sam?  159.85GB is fine, 170.81GB failed.
# However, this was before splitting out bismark2bedGraph and more aggresively removing files.
bam_size=`ls -go ${target_root}.bam | awk '{print $3}'`
# TODO: find real limit. if [ "$uncompress_bam" == "true" ] && [ $bam_size -le 171637122454 ]; then
if [ "$uncompress_bam" == "true" ] && [ $bam_size -le 200000000000 ]; then
    alignments_file=${target_root}.sam
    echo "-- Decompressing biorep bam (size: $bam_size)..."
    set -x
    samtools view ${target_root}.bam > ${alignments_file}
    set +x
    if [ $scorched == "earth" ]; then
        set -x
        rm ${target_root}.bam
        set +x
    fi
    ncores=$ncpus
else
    ncores=`expr $ncpus / 2`
    if [ "$uncompress_bam" != "true" ]; then
        echo "-- Using compressed biorep bam (size: $bam_size) and $ncores cores..."
    else
        echo "-- Using compressed bam and $ncores cores because bam_size: $bam_size exceeds limit."
    fi
fi

echo "-- Analyse methylation in ${alignments_file} and using $ncores threads..."
# NOTE: reading a bam and outputting .gz will triple the number of cores used on multi-core.
set -x
mkdir -p output/
bismark_methylation_extractor --multicore $ncores --paired-end -s --comprehensive -output output/ ${alignments_file}
mv output/*_splitting_report.txt .
mv output/${target_root}.M-bias.txt ${target_root}_mbias_report.txt
set +x

if [ $scorched == "earth" ]; then
    # Storage can be maximized by aggressively splitting out bismark2bedGraph and coverage2cytosine... and aggressively
    # removing no longer needed files.
    echo "-- Ensure cleanup to maximize available storage..."
    df -k .
    set -x
    rm -f ${alignments_file} # STORAGE IS LIMITED
    mv input/*.fa .
    rm -rf input   # Removing indexes before next step may be helpful!
    mkdir -p input
    mv *.fa input/
    ls -l output/ # Not expected to gain anything
    mv output/CpG_context_${target_root}.txt .
    mv output/CHG_context_${target_root}.txt .
    mv output/CHH_context_${target_root}.txt .
    rm -rf output
    mkdir -p output/
    mv CpG_context_${target_root}.txt output/
    mv CHG_context_${target_root}.txt output/
    mv CHH_context_${target_root}.txt output/
    set +x
    df -k .
fi

echo "-- Bismark to bedGraph..."
set -x
bismark2bedGraph --CX_context --ample_mem --dir output/ -output ${target_root}.bedGraph \
                 CpG_context_${target_root}.txt CHG_context_${target_root}.txt CHH_context_${target_root}.txt
set +x
if [ -f output/${target_root}.bedGraph ]; then
    set -x
    pigz output/${target_root}.bedGraph
    set +x
fi
set -x
mv output/${target_root}.bedGraph.gz .
set +x
if [ $scorched == "earth" ]; then
    echo "-- More scorched earth cleanup..."
    df -k .
    ls -l output/
    set -x
    mv output/${target_root}.bismark.cov.gz .
    rm -rf output
    mkdir -p output
    mv ${target_root}.bismark.cov.gz output/
    set +x
    df -k .
fi
    
echo "-- Coverage to cytosine..."
set -x
coverage2cytosine --output ${target_root}.CX_report.txt --dir 'output/' --genome 'input/' --parent_dir '/home/dnanexus' \
                  --zero --CX_context output/${target_root}.bismark.cov.gz
pigz output/${target_root}.CX_report.txt
mv output/${target_root}.CX_report.txt.gz .
set +x
if [ $scorched == "earth" ]; then
    echo "-- Final scorched earth cleanup..."
    df -k .
    ls -l output/
    rm -rf output
fi    

echo "-- The results..."
ls -l ${target_root}*
df -k .

