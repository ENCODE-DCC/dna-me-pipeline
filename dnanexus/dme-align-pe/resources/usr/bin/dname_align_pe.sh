#!/bin/bash -e

if [ $# -lt 7 ] || [ $# -gt 8 ]; then
    echo "usage v1: dname_align_pe.sh <index.tgz> <read1.fq> <read2.fq> <min_insert> <max_insert> <ncpus> <bam_root> [\"no_stats\"]"
    echo "Align pe reads with bismark.  Bowtie version is determined by index.tgz content. Is independent of DX and encodeD."
    echo "Requires cutadapt and trim_galore, bismark, bowtie1/bowtie2, and samtools on path."
    exit -1; 
fi
index_tgz=$1  # Index archive containing input/*.fa reference, input/Bisulfite_Genome/* index made with chosen bowtie version.
              # Also input/lambda/*, input/lambda/Bisulfite_Genome/*.  Will be expanded and REMOVED for storage efficiency.
read1_fq=$2   # fastq of of paired-end read1, which will be trimmed resulting in "read1_trimmed.fq.gz"
read2_fq=$3   # fastq of of paired-end read2, which will be trimmed resulting in "read2_trimmed.fq.gz"
min_insert=$4 # Minimum insert parameter (bismark -I nnn)
max_insert=$5 # Maximum insert parameter (bismark -X nnn)
ncpus=$6      # Number of cpus available (bismark uses 2 per multicore so --multi ncpus/2)
bam_root="${7}_bismark_pe"   # root name for output bam (e.g. "out" will create "out_bismark_pe.bam" and "out_bismark_pe_flagstat.txt") 

echo "-- Expect to create '${bam_root}.bam'"

echo "-- Uncompress index archive..."
set -x
tar zxvf $index_tgz
rm $index_tgz
set +x
bowtie_ver="bowtie1"
if [ -f input/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2 ]; then
    bowtie_ver="bowtie2"
fi

read1_root=${read1_fq%.gz}
read1_root=${read1_root%.fastq}
read1_root=${read1_root%.fq}
read2_root=${read2_fq%.gz}
read2_root=${read2_root%.fastq}
read2_root=${read2_root%.fq}

echo "-- Trimming reads..."
set -x
mkdir -p output/
trim_galore -o output --dont_gzip --paired $read1_fq $read2_fq
mv output/${read1_root}*.fq ${read1_root}_trimmed.fq
mv output/${read2_root}*.fq ${read2_root}_trimmed.fq
set +x
#rm -rf output/

# Note --bowtie2 and -p $nthreads are both SLOWER than single threaded bowtie1
if [ $ncpus -gt 1 ]; then
    ncores=`expr $ncpus / 2`
else
    ncores=$ncpus
fi
echo "-- Mapping to reference with bismark/${bowtie_ver} on $ncores cores..."
mkdir -p output/
if [ "$bowtie_ver" == "bowtie2" ]; then
    set -x
    bismark --bowtie2 -N 1 -L 28 --output_dir output --temp_dir output --multi $ncores \
            input -I $min_insert -X $max_insert -1 ${read1_root}_trimmed.fq -2 ${read2_root}_trimmed.fq
    set +x
else
    set -x
    bismark --bowtie1 -n 1 -l 28 --output_dir output --temp_dir output --multi $ncores \
            input -I $min_insert -X $max_insert -1 ${read1_root}_trimmed.fq -2 ${read2_root}_trimmed.fq
    set +x
fi
set -x
ls -l output/
mv output/*.bam ${bam_root}.bam
set +x

echo "-- Mapping to lambda with bismark/${bowtie_ver} on $ncores cores..."
mkdir -p output/lambda/
if [ "$bowtie_ver" == "bowtie2" ]; then
    set -x
    bismark --bowtie2 -N 1 -L 28 --output_dir output/lambda/ --temp_dir output/lambda/ --multi $ncores \
            input/lambda -I $min_insert -X $max_insert -1 ${read1_root}_trimmed.fq -2 ${read2_root}_trimmed.fq
    set +x
else
    set -x
    mkdir -p output/lambda/
    bismark --bowtie1 -n 1 -l 28 --output_dir output/lambda/ --temp_dir output/lambda/ --multi $ncores \
            input/lambda -I $min_insert -X $max_insert -1 ${read1_root}_trimmed.fq -2 ${read2_root}_trimmed.fq
    set +x
fi
# Don't care about the bam, only the map_report

echo "-- Combine map reports..."
echo "===== bismark reference =====" > ${bam_root}_map_report.txt
cat output/*PE_report.txt           >> ${bam_root}_map_report.txt
echo " "                            >> ${bam_root}_map_report.txt
echo "===== bismark lambda ====="   >> ${bam_root}_map_report.txt
cat output/lambda/*PE_report.txt    >> ${bam_root}_map_report.txt

if [ $# -lt 8 ]; then
    echo "-- Collect bam stats..."
    set -x
    samtools flagstat ${bam_root}.bam > ${bam_root}_flagstat.txt
    samtools stats ${bam_root}.bam > ${bam_root}_samstats.txt
    head -3 ${bam_root}_samstats.txt
    grep ^SN ${bam_root}_samstats.txt | cut -f 2- > ${bam_root}_samstats_summary.txt
    set +x
fi

echo "-- The results..."
ls -l ${target_root}*
df -k .

