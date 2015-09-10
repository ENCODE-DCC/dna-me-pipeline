#!/bin/bash
# dme-map-se.sh - WGBS ENCODE Pipeline step: Trim and map single ended reads for technical replicate using Bismark

set -x
set +e

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of reads:     '${reads[@]}'"
    echo "* Value of dme_ix:    '$dme_ix'"
    #echo "* Value of min_insert:     $min_insert"
    #echo "* Value of max_insert:     $max_insert"

    # NOTE: not expecting an array of files, but supporting it nonetheless.
    # load and concat reads
    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads[@]}
    do
        file_root=`dx describe "${reads[$ix]}" --name`
        file_root=${file_root%.fastq.gz}
        file_root=${file_root%.fq.gz}
        if [ "${outfile_name}" == "" ]; then
            outfile_name="${file_root}"
        else
            outfile_name="${file_root}_${outfile_name}"
            if [ "${concat}" == "" ]; then
                outfile_name="${outfile_name}_concat" 
                concat="s concatenated as"
            fi
        fi
        echo "* Downloading and concatenating ${file_root}.fq.gz file..."
        dx download "${reads[$ix]}" -o - | gunzip >> concat.fq
    done
    # Try to simplify the names
    rep_root=""
    if [ -f /usr/bin/parse_property.py ]; then
        rep_root=`parse_property.py -f "'${reads1[0]}'" --project "${DX_PROJECT_CONTEXT_ID}" --root_name --quiet`
    fi
    if [ "$rep_root" != "" ]; then
        outfile_name="${rep_root}_reads"
    else
        outfile_name="${outfile_name}_reads"
    fi
    mv concat.fq ${outfile_name}.fq
    #echo "* Gzipping file..."
    #gzip ${outfile_name}.fq
    reads_root=${outfile_name}
    echo "* Fastq${concat} file: '${reads_root}.fq'"
    ls -l ${reads_root}.fq

    echo "* Download and uncompress index archive..."
    dx download "$dme_ix" -o - | tar zxvf -
    #ls -l input

    bam_root="${reads_root}_bismark_techrep"
    # Try to simplify the names
    if [ "$rep_root" != "" ]; then
        bam_root="${rep_root}_bismark_techrep"
    fi
    echo "* Expect to create '${bam_root}.bam'"

    echo "* Trimming reads..."
    # NOTE: 2 different versions of mit-trim.py !  
    set -x
    mott-trim-se.py -q 3 -m 30 -t sanger ${reads_root}.fq > ${reads_root}_trimmed.fq
    set +x

    echo "* Mapping to reference with bismark/bowtie..."
    set -x
    mkdir output
    bismark -n 1 -l 28 -output_dir output --temp_dir output input ${reads_root}_trimmed.fq
    set +x

    echo "* Compressing with samtools..."
    set -x
    samtools view -Sb -@ ${nthreads} output/${reads_root}_trimmed.fq_bismark_se.sam > ${bam_root}.bam
    cat output/*SE_report.txt > ${bam_root}_ref_report.txt
    set +x

    echo "* Mapping to lambda with bismark/bowtie..."
    set -x
    mkdir output
    bismark -n 1 -l 28 -output_dir output/lambda --temp_dir output/lambda input/lambda ${reads_root}_trimmed.fq
    samtools view -Sb -@ ${nthreads} output/lambda/${reads_root}_trimmed.fq_bismark_se.sam > ${bam_root}_lambda.bam
    cat output/lambda/*SE_report.txt > ${bam_root}_lambda_report.txt
    set +x

    echo "* Collect bam stats..."
    set -x
    samtools flagstat ${bam_root}.bam > ${bam_root}_flagstat.txt
    samtools stats ${bam_root}.bam > ${bam_root}_samstats.txt
    head -3 ${bam_root}_samstats.txt
    grep ^SN ${bam_root}_samstats.txt | cut -f 2- > ${bam_root}_samstats_summary.txt
    set +x

    echo "* Prepare metadata..."
    echo "===== bismark reference =====" > ${bam_root}_map_report.txt
    cat ${bam_root}_ref_report.txt      >> ${bam_root}_map_report.txt
    echo " "                            >> ${bam_root}_map_report.txt
    echo "===== bismark lambda ====="   >> ${bam_root}_map_report.txt
    cat ${bam_root}_lambda_report.txt   >> ${bam_root}_map_report.txt

    qc_stats=''
    reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n bismark_ref -f ${bam_root}_map_report.txt`
        meta=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_flagstat.txt`
        qc_stats=`echo $qc_stats, $meta`
        reads=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_flagstat.txt -k total`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_root}_samstats_summary.txt`
        read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_root}_samstats_summary.txt -k "average length"`
        qc_stats=`echo $qc_stats, $meta`
    fi
    # All qc to one file per target file:
    echo "===== samtools flagstat =====" > ${bam_root}_qc.txt
    cat ${bam_root}_flagstat.txt        >> ${bam_root}_qc.txt
    echo " "                            >> ${bam_root}_qc.txt
    echo "===== samtools stats ====="   >> ${bam_root}_qc.txt
    cat ${bam_root}_samstats.txt        >> ${bam_root}_qc.txt

    echo "* Upload results..."
    ls -l /home/dnanexus/output
    #echo `ls /home/dnanexus/output`
    #outfile="$read_fn".mapped_methylseq.tgz
    #tar zcvf $outfile output
    #mapped_files=$(dx upload /home/dnanexus/$outfile --brief)
    bam_techrep=$(dx upload /home/dnanexus/${bam_root}.bam --details "{ $qc_stats }" --property SW="$versions" \
                                               --property reads="$reads" --property read_length="$read_len" --brief)
    bam_techrep_qc=$(dx upload ${bam_root}_qc.txt --details "{ $qc_stats }" --property SW="$versions" --brief)
    map_techrep=$(dx upload ${bam_root}_map_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_techrep "$bam_techrep" --class=file
    dx-jobutil-add-output bam_techrep_qc "$bam_techrep_qc" --class=file
    dx-jobutil-add-output map_techrep "$map_techrep" --class=file

    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
 }
