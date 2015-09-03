#!/bin/bash
# dme-extract.sh - WGBS ENCODE Pipeline step: Extract methylation and report Whole Genome Bisulphite Analysis

set -x
set +e

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo " Value of gzip: '$gzip'"

    echo "getting files"
    dx download "$genome" -o - | gunzip > genome.fa
    mapped_fn=`dx describe "$mapped_files" --name | cut -d'.' -f1`
    dx download "$mapped_files" -o - | tar zxvf -

    dx download "$chrom_sizes" -o chrom.sizes

    mkdir input
    mv genome.fa input
    echo "Analyse methylation"
    outfile="$mapped_fn".fq_bismark

    gzipFlag=""
    if [ "$gzip" == "true" ]; then
        echo '* Adding gzip flag'
        gzipFlag="--gzip"
    else
        echo '* No gzip flaf'
    fi


    bismark_methylation_extractor "$gzipFlag" -s --comprehensive --cytosine_report --CX_context --ample_mem\
      --output /home/dnanexus/output/ --zero_based --genome_folder input output/"$outfile".sam

    samtools view -Sb output/"$outfile".sam > output/"$outfile".bam
    echo "Creat QC reports"
    cxrepo-bed.py -o /home/dnanexus/output /home/dnanexus/output/"$mapped_fn".fq_bismark.CX_report.txt

    # Fill in your application code here.
    #
    # To report any recognized errors in the correct format in
    # $HOME/job_error.json and exit this script, you can use the
    # dx-jobutil-report-error utility as follows:
    #
    #   dx-jobutil-report-error "My error message"
    #
    # Note however that this entire bash script is executed with -e
    # when running in the cloud, so any line which returns a nonzero
    # exit code will prematurely exit the script; if no error was
    # reported in the job_error.json file, then the failure reason
    # will be AppInternalError with a generic error message.

    # The following line(s) use the dx command-line tool to upload your file
    # outputs after you have created them on the local file system.  It assumes
    # that you have used the output field name for the filename for each output,
    # but you can change that behavior to suit your needs.  Run "dx upload -h"
    # to see more options to set metadata.

    echo `ls /home/dnanexus/output`
    mv /home/dnanexus/output/CG_"$mapped_fn".fq_bismark.CX_report.txt "$mapped_fn"_CG_bismark.bed
    mv /home/dnanexus/output/CHG_"$mapped_fn".fq_bismark.CX_report.txt "$mapped_fn"_CHG_bismark.bed
    mv /home/dnanexus/output/CHH_"$mapped_fn".fq_bismark.CX_report.txt "$mapped_fn"_CHH_bismark.bed

    echo "Convert to BigBed"
    bedToBigBed "$mapped_fn"_CG_bismark.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes "$mapped_fn"_CG_bismark.bb
    bedToBigBed "$mapped_fn"_CHG_bismark.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes "$mapped_fn"_CHG_bismark.bb
    bedToBigBed "$mapped_fn"_CHH_bismark.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes "$mapped_fn"_CHH_bismark.bb

    gzip *.bed
    echo "Uploading files"
    find
    CG=$(dx upload "$mapped_fn"_CG_bismark.bed.gz --brief)
    CHG=$(dx upload "$mapped_fn"_CHG_bismark.bed.gz --brief)
    CHH=$(dx upload "$mapped_fn"_CHH_bismark.bed.gz --brief)

    CGbb=$(dx upload "$mapped_fn"_CG_bismark.bb --brief)
    CHGbb=$(dx upload "$mapped_fn"_CHG_bismark.bb --brief)
    CHHbb=$(dx upload "$mapped_fn"_CHH_bismark.bb --brief)

    mapped_reads=$(dx upload /home/dnanexus/output/"$outfile".bam --brief)
    cat output/*E_report.txt > output/$mapped_fn.fq_bismark_map_report.txt
    map_report=$(dx upload /home/dnanexus/output/"$mapped_fn".fq_bismark_map_report.txt --brief)
    M_bias_report=$(dx upload /home/dnanexus/output/"$mapped_fn".fq_bismark.M-bias.txt --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.
    echo "Adding output -- files should be renamed"

    dx-jobutil-add-output CG "$CG" --class=file
    dx-jobutil-add-output CHG "$CHG" --class=file
    dx-jobutil-add-output CHH "$CHH" --class=file
    dx-jobutil-add-output CGbb "$CGbb" --class=file
    dx-jobutil-add-output CHGbb "$CHGbb" --class=file
    dx-jobutil-add-output CHHbb "$CHHbb" --class=file
    dx-jobutil-add-output mapped_reads "$mapped_reads" --class=file
    dx-jobutil-add-output map_report "$map_report" --class=file
    dx-jobutil-add-output M_bias_report "$M_bias_report" --class=file
}
