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

    echo "* Download and untar index file..."
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

    echo "* Mapping with bismark/bowtie..."
    set -x
    mkdir output
    bismark -n 1 -l 28 -output_dir output --temp_dir output input ${reads_root}_trimmed.fq
    set +x

    echo "* Compressing with samtools..."
    set -x
    samtools -Sb -@ ${nthreads} output/${reads_root}.fq_bismark_se.sam ${bam_root}.bam
    set +x

    echo "* Prepare metadata..."
    qc_stats=''
    reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_flagstat.txt`
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

    dx-jobutil-add-output bam_techrep "$bam_techrep" --class=file
    dx-jobutil-add-output bam_techrep_qc "$bam_techrep_qc" --class=file

    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
 }
