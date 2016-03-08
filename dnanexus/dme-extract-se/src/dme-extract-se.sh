#!/bin/bash
# dme-extract-se.sh - WGBS ENCODE Pipeline step: Merge techrep bams, extract single-end methylation and report WGBS Analysis.

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_set:        '$bam_set'"
    echo "* Value of map_report_set: '$map_report_set'"
    echo "* Value of dme_ix:         '$dme_ix'"
    echo "* Value of uncompress_bam: '$uncompress_bam'"

    # NOTE: dme-align produces *_techrep_bismark.bam and dme-extract merges 1+ techrep bams into a *_bismark_biorep.bam.
    #       The reason for the name 'word' order is so thal older *_bismark.bam alignments are recognizable as techrep bams

    target_root=""
    merged=""
    tech_reps=""
    for ix in ${!bam_set[@]}
    do
        file_root=`dx describe "${bam_set[$ix]}" --name`
        file_root=${file_root%_techrep_bismark.bam}
        file_root=${file_root%_bismark.bam}
        if [ "${target_root}" == "" ]; then
            target_root="${file_root}"
        else
            target_root="${file_root}_${target_root}"
            if [ "${merged}" == "" ]; then
                target_root="${target_root}_bismark_biorep" 
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
        echo "* Downloading ${file_root}_techrep_bismark.bam file..."
        dx download "${bam_set[$ix]}" -o ${file_root}_techrep_bismark.bam
        if [ ! -e sofar.bam ]; then
            mv ${file_root}_techrep_bismark.bam sofar.bam
        else
            echo "* Merging..."
            # NOTE: keeps the first header
            set -x
            samtools cat sofar.bam ${file_root}_techrep_bismark.bam > merging.bam
            mv merging.bam sofar.bam
            rm ${file_root}_techrep_bismark.bam # STORAGE IS LIMITED
            set +x
        fi
    done
    if [ "$exp_id" != "" ] && [ "$tech_reps" != "" ]; then
        target_root="${exp_id}_${tech_reps}_bismark_biorep"
    fi
    echo "* Merged alignments file will be: '${target_root}.bam'"

    # At this point there is a 'sofar.bam' with one or more input bams
    if [ "${merged}" == "" ]; then
        target_root="${file_root}_bismark_biorep"
        set -x
        mv sofar.bam ${target_root}.bam
        set +x
        echo "* Only one input file, no merging required."
    else
        # sorting needed due to samtools cat
        echo "* Sorting merged bam..."
        set -x
        samtools sort -@ 16 2900M -f sofar.bam sorted.bam
        mv sorted.bam ${target_root}.bam
        rm sofar.bam # STORAGE IS LIMITED
        set +x
        echo "* Files merged into '${target_root}.bam'"
    fi 

    # Working on map_reports now
    all_reports=""
    echo "### Combined Bismark map report for several technical replicates ###" > ${target_root}_map_report.txt
    echo " " >> ${target_root}_map_report.txt
    for ix in ${!map_report_set[@]}
    do
        file_root=`dx describe "${map_report_set[$ix]}" --name`
        file_root=${file_root%_techrep_bismark_map_report.txt}
        file_root=${file_root%_bismark_map_report.txt}
        file_root=${file_root%_map_report.txt}
        echo "###################################" >> ${target_root}_map_report.txt
        echo "### Map report for ${file_root} ###" >> ${target_root}_map_report.txt
        echo "* Downloading ${file_root}_techrep_bismark_map_report.txt file..."
        dx download "${map_report_set[$ix]}" -o ${file_root}_map_report.txt
        cat ${file_root}_map_report.txt >> ${target_root}_map_report.txt
        if [ "${all_reports}" == "" ]; then
            all_reports="${file_root}_map_report.txt"
        else
            all_reports="${all_reports},${file_root}_map_report.txt"
        fi
    done
    if [ "${all_reports}" == "${file_root}_map_report.txt" ]; then # only one
        cp ${file_root}_map_report.txt ${target_root}_map_report.txt
    fi
    
    echo "* Download index archive..."
    dx download "$dme_ix" -o index.tgz

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dname_extract_se.sh index.tgz ${target_root}.bam 32 $uncompress_bam --scorched_earth
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="

    echo "* Prepare metadata..."
    qc_stats=''
    reads=0
    #read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n bismark_map -f ${all_reports}`
        meta=`qc_metrics.py -n samtools_flagstats -f ${target_root}_flagstat.txt`
        qc_stats=`echo $qc_stats, $meta`
        reads=`qc_metrics.py -n samtools_flagstats -f ${target_root}_flagstat.txt -k total`
        meta=`qc_metrics.py -n bismark_extract -f *_splitting_report.txt`
        qc_stats=`echo $qc_stats, $meta`
    fi
    # All qc to one file per target file:
    echo "===== samtools flagstat ====="  > ${target_root}_qc.txt
    cat ${target_root}_flagstat.txt      >> ${target_root}_qc.txt
    echo " "                                                           >> ${target_root}_qc.txt
    echo "===== bismark_methylation_extractor: splitting_report =====" >> ${target_root}_qc.txt
    cat *_splitting_report.txt                                         >> ${target_root}_qc.txt

    echo "* Uploading files..."
    ls -l 
    # NOTE: Not saving merged bam
    #bam_biorep=$(dx upload ${target_root}.bam --details "{ $qc_stats }" --property SW="$versions" \
    #                                           --property reads="$reads" --property read_length="$read_len" --brief)
    bam_biorep_qc=$(dx upload ${target_root}_qc.txt --details "{ $qc_stats }" --property SW="$versions" --brief)
    map_biorep=$(dx upload ${target_root}_map_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)
    mbias_report=$(dx upload ${target_root}_mbias_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)

    #dx-jobutil-add-output bam_biorep "$bam_biorep" --class=file
    dx-jobutil-add-output bam_biorep_qc "$bam_biorep_qc" --class=file
    dx-jobutil-add-output map_biorep "$map_biorep" --class=file
    dx-jobutil-add-output mbias_report "$mbias_report" --class=file

    # Interim files
    cx_report=$(dx upload ${target_root}.CX_report.txt.gz --details "{ $qc_stats }" --property SW="$versions" --brief)
    bg_gz=$(dx upload ${target_root}.bedGraph.gz --details "{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output cx_report "$cx_report" --class=file
    dx-jobutil-add-output bg_gz "$bg_gz" --class=file

    # Metadata
    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
