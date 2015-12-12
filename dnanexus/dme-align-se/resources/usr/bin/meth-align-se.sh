#!/bin/bash -e

if [ $# -ne 4 ]; then
    echo "usage v1: meth-align-se.sh <index.tgz> <reads.fq> <ncpus> <bam_root>"
    echo "Align se reads with bismark.  Bowtie version is determined by index.tgz content. Is independent of DX and encodeD."
    echo "Requires cutadapt and trim_galore, bismark, bowtie1/bowtie2, and samtools on path."
    exit -1; 
fi
index_tgz=$1  # Index archive containing input/*.fa reference, input/Bisulfite_Genome/* index made with chosen bowtie version.
              # Also input/lambda/*, input/lambda/Bisulfite_Genome/*.  Will be expanded and REMOVED for storage efficiency.
reads_fq=$2   # fastq of of paired-end read1, which will be trimmed resulting in "read1_trimmed.fq.gz"
ncpus=$3      # Number of cpus available (bismark uses 2 per multicore so --multi ncpus/2)
bam_root=$4   # root name for output bam (e.g. "out_bam" will create "out_bam_pe.bam" and "out_bam_pe_flagstats.txt") 

echo "-- Uncompress index archive..."
set -x
tar zxvf $index_tgz
rm $index_tgz
set +x
bowtie_ver="bowtie1"
if [ -f input/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2 ]; then
    bowtie_ver="bowtie2"
fi

reads_root=${reads_fq%.gz}
reads_root=${reads_root%.fastq}
reads_root=${reads_root%.fq}

echo "-- Trimming reads..."
set -x
mkdir -p output/
trim_galore -o output --dont_gzip $reads_fq
mv output/${reads_root}*.fq ${reads_root}_trimmed.fq
set +x
#rm -rf output/
###echo "-- Normalize fastq name and uncompressed if necessary..."
###if [[ ${reads_fq} =~ \.gz$ ]]; then
###    set -x
###    if [ ${reads_fq} != ${reads_root}.fq.gz ]; then
###        mv ${reads_fq} ${reads_root}.fq.gz
###    fi
###    gunzip ${reads_root}.fq.gz
###    set +x
###elif [ ${reads_fq} != ${reads_root}.fq ]; then
###    set -x
###    mv ${reads_fq} ${reads_root}.fq
###    set +x
###fi
###echo "-- Trimming reads..."
###set -x
###mott-trim-se.py -q 3 -m 30 -t sanger ${reads_root}.fq > ${reads_root}_trimmed.fq
###set +x

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
            input ${reads_root}_trimmed.fq
    set +x
else
    set -x
    bismark --bowtie1 -n 1 -l 28 --output_dir output --temp_dir output --multi $ncores \
            input ${reads_root}_trimmed.fq
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
            input/lambda ${reads_root}_trimmed.fq
    set +x
else
    set -x
    mkdir -p output/lambda/
    bismark --bowtie1 -n 1 -l 28 --output_dir output/lambda/ --temp_dir output/lambda/ --multi $ncores \
            input/lambda ${reads_root}_trimmed.fq
    set +x
fi
# Don't care about the bam, only the map_report

echo "-- Combine map reports..."
echo "===== bismark reference =====" > ${bam_root}_map_report.txt
cat output/*SE_report.txt           >> ${bam_root}_map_report.txt
echo " "                            >> ${bam_root}_map_report.txt
echo "===== bismark lambda ====="   >> ${bam_root}_map_report.txt
cat output/lambda/*SE_report.txt    >> ${bam_root}_map_report.txt

echo "-- Collect bam stats..."
set -x
samtools flagstat ${bam_root}.bam > ${bam_root}_flagstat.txt
samtools stats ${bam_root}.bam > ${bam_root}_samstats.txt
head -3 ${bam_root}_samstats.txt
grep ^SN ${bam_root}_samstats.txt | cut -f 2- > ${bam_root}_samstats_summary.txt
set +x

echo "-- The results..."
ls -l ${target_root}*
df -k .

