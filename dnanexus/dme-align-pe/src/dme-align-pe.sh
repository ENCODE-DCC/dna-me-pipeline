#!/bin/bash
# dme-align-pe.sh - align pe reads with bismark/bowtie

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        index_file=`dx describe "$dme_ix" --name`
        if [[ $index_file =~ *"bowtie1"* ]]; then
            versions=`tool_versions.py --applet dme-align-pe-bowtie1`
        else
            versions=`tool_versions.py --dxjson dnanexus-executable.json`
        fi
    fi

    echo "* Value of pair 1 reads: '${reads1[@]}"
    echo "* Value of pair 2 reads: '${reads2[@]}"
    echo "* Value of dme_ix:       '$dme_ix'"
    echo "* Value of ncpus:         $ncpus"
    echo "* Value of min_insert:    $min_insert"
    echo "* Value of max_insert:    $max_insert"

    # NOTE: dme-align produces *_techrep_bismark.bam and dme-extract merges 1+ techrep bams into a *_bismark_biorep.bam.
    #       The reason for the name 'word' order is so thal older *_bismark.bam alignments are recognizable as techrep bams

    # NOTE: not expecting an array of files, but supporting it nonetheless.
    # load and concat reads1
    exp_rep_root=""
    if [ -f /usr/bin/parse_property.py ]; then
        new_root=`parse_property.py -f "'${reads1[0]}'" --project "${DX_PROJECT_CONTEXT_ID}" --root_name --quiet`
        if [ "$new_root" != "" ]; then
            exp_rep_root="${new_root}"
        fi
    fi
    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads1[@]}
    do
        file_root=`dx describe "${reads1[$ix]}" --name`
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
        echo "* Downloading concatenating ${file_root}.fq.gz file..."
        dx download "${reads1[$ix]}" -o - | gunzip >> concat.fq
    done
    # Try to simplify the names
    if [ "${concat}" != "" ]; then
        if [ "${exp_rep_root}" != "" ]; then
            outfile_name="${exp_rep_root}_reads1"
        elif [ ${#outfile_name} -gt 200 ]; then
            outfile_name="concatenated_reads1"
        fi
    fi
    mv concat.fq ${outfile_name}.fq
    #echo "* Gzipping file..."
    #gzip ${outfile_name}.fq
    reads1_root=${outfile_name}
    echo "* Reads1 fastq${concat} file: '${reads1_root}.fq'"
    ls -l ${reads1_root}.fq

    # load and concat reads2
    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads2[@]}
    do
        file_root=`dx describe "${reads2[$ix]}" --name`
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
        dx download "${reads2[$ix]}" -o - | gunzip >> concat.fq
    done
    # Try to simplify the names
    if [ "${concat}" != "" ]; then
        if [ "${exp_rep_root}" != "" ]; then
            outfile_name="${exp_rep_root}_reads2"
        elif [ ${#outfile_name} -gt 200 ]; then
            outfile_name="concatenated_reads2"
        fi
    fi
    mv concat.fq ${outfile_name}.fq
    #echo "* Gzipping file..."
    #gzip ${outfile_name}.fq
    reads2_root=${outfile_name}
    echo "* Reads2 fastq${concat} file: '${reads2_root}.fq'"
    ls -l ${reads2_root}.fq

    echo "* Download index archive..."
    dx download "$dme_ix" -o index.tgz

    bam_root="${reads1_root}_${reads2_root}_techrep"
    # Try to simplify the names
    if [ "$exp_rep_root" != "" ]; then
        bam_root="${exp_rep_root}_techrep"
    fi

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dname_align_pe.sh index.tgz ${reads1_root}.fq ${reads2_root}.fq $min_insert $max_insert $ncpus $bam_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    bam_root="${bam_root}_bismark_pe"

    echo "* Prepare metadata..."
    qc_stats=''
    reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n bismark_map -f ${bam_root}_map_report.txt`
        meta=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_flagstat.txt`
        qc_stats=`echo $qc_stats, $meta`
        reads=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_flagstat.txt -k total`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_root}_samstats_summary.txt`
        read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_root}_samstats_summary.txt -k "average length"`
        qc_stats=`echo $qc_stats, $meta`
    fi
    # All qc to one file:
    cat ${bam_root}_map_report.txt      >> ${bam_root}_qc.txt
    echo " "                            >> ${bam_root}_qc.txt
    echo "===== samtools flagstat =====" > ${bam_root}_qc.txt
    cat ${bam_root}_flagstat.txt        >> ${bam_root}_qc.txt
    echo " "                            >> ${bam_root}_qc.txt
    echo "===== samtools stats ====="   >> ${bam_root}_qc.txt
    cat ${bam_root}_samstats.txt        >> ${bam_root}_qc.txt

    echo "* Upload results..."
    #ls -l /home/dnanexus/output
    bam_techrep=$(dx upload ${bam_root}.bam --details "{ $qc_stats }" --property SW="$versions" \
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
