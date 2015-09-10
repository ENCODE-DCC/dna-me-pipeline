#!/bin/bash
# dme-merge-bams.sh - WGBS ENCODE Pipeline step: Merge two or more technical replicate bams for the ENCODE DNase-seq pipeline

main() {
    # Executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_set:        '$bam_set'"
    echo "* Value of map_report_set: '$bam_set'"
    echo "* Value of nthreads:       '$nthreads'"
    
    outfile_name=""
    merged=""
    tech_reps=""
    for ix in ${!bam_set[@]}
    do
        file_root=`dx describe "${bam_set[$ix]}" --name`
        file_root=${file_root%_bismark_techrep.bam}
        file_root=${file_root%_bismark.bam}
        if [ "${outfile_name}" == "" ]; then
            outfile_name="${file_root}"
        else
            outfile_name="${file_root}_${outfile_name}"
            if [ "${merged}" == "" ]; then
                outfile_name="${outfile_name}_bismark_biorep" 
                merged="s merged as"
            fi
        fi
        # Try to simplify the names
        if [ -f /usr/bin/parse_property.py ]; then
            if [ "$exp_id" == "" ]; then
                exp_id=`parse_property.py -f "'${bam_set[$ix]}'" --project "${DX_PROJECT_CONTEXT_ID}" --exp_id`
            fi
            rep_tech=`parse_property.py -f "'${bam_set[$ix]}'" --project "${DX_PROJECT_CONTEXT_ID}" --rep_tech`
            if [ "$rep_tech" != "" ]; then
                if  [ "$tech_reps" != "" ]; then
                    tech_reps="${tech_reps}_${rep_tech}"
                else
                    tech_reps="${rep_tech}"
                fi
            fi
        fi
        echo "* Downloading ${file_root}_bismark_techrep.bam file..."
        dx download "${bam_set[$ix]}" -o ${file_root}_bismark_techrep.bam
        if [ ! -e sofar.bam ]; then
            mv ${file_root}_bismark_techrep.bam sofar.bam
        else
            echo "* Merging..."
            # NOTE: keeps the first header
            set -x
            samtools cat sofar.bam ${file_root}_bismark_techrep.bam > merging.bam
            mv merging.bam sofar.bam
            set +x
        fi
    done
    if [ "$exp_id" != "" ] && [ "$tech_reps" != "" ]; then
        outfile_name="${exp_id}_${tech_reps}_bismark_biorep"
    fi
    echo "* Merged alignments file will be: '${outfile_name}.bam'"

    # Working on map_reports now
    all_reports=""
    echo "### Combined Bismark map report for several technical replicates ###" > ${outfile_name}_map_report.txt
    echo " " >> ${outfile_name}_map_report.txt
    for ix in ${!map_report_set[@]}
    do
        file_root=`dx describe "${map_report_set[$ix]}" --name`
        file_root=${file_root%_bismark_techrep_map_report.txt}
        file_root=${file_root%_map_report.txt}
        echo "###################################" >> ${outfile_name}_map_report.txt
        echo "### Map report for ${file_root} ###" >> ${outfile_name}_map_report.txt
        echo "* Downloading ${file_root}_bismark_techrep_map_report.txt file..."
        dx download "${bam_set[$ix]}" -o ${file_root}_map_report.txt
        cat ${file_root}_map_report.txt >> ${outfile_name}_map_report.txt
        if [ "${all_reports}" == "" ]; then
            all_reports="${file_root}_map_report.txt"
        else
            all_reports="${all_reports},${file_root}_map_report.txt"
        fi
    done
    if [ "${all_reports}" == "${file_root}_map_report.txt" ]; then # only one
        cp ${file_root}_map_report.txt ${outfile_name}_map_report.txt
    fi
    qc_stats=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n bismark_map -f ${all_reports}`
    fi
    
    # TODO: sorting needed?
    echo "* Sorting merged bam..."
    set -x
    samtools sort -@ $nthreads -m 6G -f sofar.bam sorted
    samtools view -hb sorted.bam > sofar.bam 
    set +x
    
    # At this point there is a 'sofar.bam' with one or more input bams
    if [ "${merged}" == "" ]; then
        outfile_name="${file_root}_bismark_biorep"
        mv sofar.bam ${outfile_name}.bam
        echo "* Only one input file, no merging required."
    else
        mv sofar.bam ${outfile_name}.bam
        echo "* Files merged into '${outfile_name}.bam'"
    fi 

    echo "* Collect bam stats..."
    set -x
    samtools flagstat ${outfile_name}.bam > ${outfile_name}_flagstat.txt
    samtools stats ${outfile_name}.bam > ${outfile_name}_samstats.txt
    head -3 ${outfile_name}_samstats.txt
    grep ^SN ${outfile_name}_samstats.txt | cut -f 2- > ${outfile_name}_samstats_summary.txt
    set +x


    echo "* Prepare metadata..."
    reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n samtools_flagstats -f ${outfile_name}_flagstat.txt`
        qc_stats=`echo $qc_stats, $meta`
        reads=`qc_metrics.py -n samtools_flagstats -f ${outfile_name}_flagstat.txt -k total`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${outfile_name}_samstats_summary.txt`
        read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${outfile_name}_samstats_summary.txt -k "average length"`
        qc_stats=`echo $qc_stats, $meta`
    fi
    # All qc to one file per target file:
    echo "===== samtools flagstat =====" > ${outfile_name}_qc.txt
    cat ${outfile_name}_flagstat.txt    >> ${outfile_name}_qc.txt
    echo " "                            >> ${outfile_name}_qc.txt
    echo "===== samtools stats ====="   >> ${outfile_name}_qc.txt
    cat ${outfile_name}_samstats.txt    >> ${outfile_name}_qc.txt

    echo "* Upload results..."
    bam_biorep=$(dx upload ${outfile_name}.bam --details "{ $qc_stats }" --property SW="$versions" \
                                               --property reads="$reads" --property read_length="$read_len" --brief)
    bam_biorep_qc=$(dx upload ${outfile_name}_qc.txt --details "{ $qc_stats }" --property SW="$versions" --brief)
    map_biorep=$(dx upload ${outfile_name}_map_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_biorep "$bam_biorep" --class=file
    dx-jobutil-add-output bam_biorep_qc "$bam_biorep_qc" --class=file
    dx-jobutil-add-output map_biorep "$map_biorep" --class=file

    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
