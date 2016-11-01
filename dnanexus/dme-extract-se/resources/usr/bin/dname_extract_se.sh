#!/bin/bash -e

if [ $# -lt 4 ] || [ $# -gt 6 ]; then
    echo "usage v1: dname_extract_se.sh <index.tgz> <bismark.bam> <ncpus> [uncompress_bam] [dedup] [--scorched_earth]"
    echo "Extracts methylation from single-end bismark bam.  Is independent of DX and ENCODE."
    echo "If --scorched_earth will remove everything, including input bam and index in order to maximize available storage."
    echo "If --dedup will run bismark PCR duplicate remover"
    exit -1; 
fi
index_tgz=$1  # Index archive containing input/*.fa reference, input/Bisulfite_Genome/* index made with chosen bowtie version.
               # Also input/lambda/*, input/lambda/Bisulfite_Genome/*.  Will be expanded and REMOVED for storage efficiency
bismark_bam=$2  # bam input aligned with bismark and bowtie version that matches index_tgz.
ncpus=$3         # Number of cores available for bismark --multi
uncompress_bam=$4 # Uncompress bam if possible: trade off: storage vs threads available during bismark
dedup=$5
scorched="no"   # --scorched_earth will remove everything, including input bam and index in order to maximize available storage.
if [ $# -eq 6 ] && [ "$6" == "--scorched_earth" ]; then
    echo "-- Scorched earth policy to maximize available storage..."
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
# TODO: find real limit. if [ $scorched != "earth" ] && [ $bam_size -le 171637122454 ]; then
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
    if [ "$uncompress_bam" == "true" ]; then
        echo "-- Using compressed bam and $ncores cores because bam_size: $bam_size exceeds limit."
    else
        echo "-- Using compressed biorep bam (size: $bam_size) and $ncores cores..."
    fi
fi

if [ $scorched == "earth" ]; then
    # Storage can be maximized by aggressively splitting out bismark2bedGraph and coverage2cytosine... and aggressively
    # removing no longer needed files.
    echo "-- Scorched earth means remove index, only keeping the *.fa file..."
    df -k .
    set -x
    mv input/*.fa .
    rm -rf input   # Removing indexes before next step may be helpful!
    mkdir -p input
    mv *.fa input/
    pigz input/*.fa
    set +x
    df -k .
fi

echo "-- Analyse methylation in ${alignments_file} and using $ncores threads..."
# NOTE: reading a bam and outputting .gz will triple the number of cores used on multi-core.
set -x
mkdir -p output/

if [ $dedup == "true" ]; then
    echo "-- Deduplicating reads"
    ### from HAIB
    ### Run the deduplication, and remove the pcr duplicates from unsorted_bam_files (i.e the sequence aligning to the same genomic positions).
    deduplicate_bismark -s --barcode ${alignments_file}
fi

bismark_methylation_extractor --multicore $ncores --single-end --comprehensive -output output/ ${alignments_file}


mv output/*_splitting_report.txt .
mv output/${target_root}.M-bias.txt ${target_root}_mbias_report.txt
set +x

txt_suffix="txt"
if [ $scorched == "earth" ]; then
    # Storage can be maximized by splitting out bismark2bedGraph and coverage2cytosine... and aggressively
    # removing no longer needed files.
    echo "-- Scorched earth cleanup to maximize available storage..."
    df -k .
    set -x
    rm -f ${alignments_file} # STORAGE IS LIMITED
    ls -l output/ # Not expected to gain anything
    mv output/CpG_context_${target_root}.txt .
    mv output/CHG_context_${target_root}.txt .
    mv output/CHH_context_${target_root}.txt .
    rm -rf output
    mkdir -p output/
    mv CpG_context_${target_root}.txt output/
    mv CHG_context_${target_root}.txt output/
    mv CHH_context_${target_root}.txt output/
    pigz output/CpG_context_${target_root}.txt
    pigz output/CHG_context_${target_root}.txt
    pigz output/CHH_context_${target_root}.txt
    txt_suffix="txt.gz"
    set +x
    df -k .
fi

echo "-- Bismark to bedGraph..."
set -x
bismark2bedGraph --CX_context --ample_mem --dir output/ -output ${target_root}.bedGraph \
                 CpG_context_${target_root}.${txt_suffix} CHG_context_${target_root}.${txt_suffix} \
                 CHH_context_${target_root}.${txt_suffix}
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
    echo "-- Must expand genome fasta file, though..."
    # Must expand input/*.fa.gz, though
    set -x
    gunzip input/*.fa.gz
    set +x    
    df -k .
fi
    
echo "-- Coverage to cytosine..."
set -x
coverage2cytosine --gzip --output ${target_root}.CX_report.txt.gz --dir 'output/' --genome_folder 'input/' --parent_dir '/home/dnanexus' \
                  --zero_based --CX_context ${target_root}.bismark.cov.gz
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

