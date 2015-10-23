#!/bin/bash
# dme-extract-pe.sh - WGBS ENCODE Pipeline step: Extract paired-end methylation and report Whole Genome Bisulphite Analysis.

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_set:        '$bam_set'"
    echo "* Value of map_report_set: '$bam_set'"
    echo "* Value of dme_ix:         '$dme_ix'"
    echo "* Value of nthreads:       '$nthreads'"
    
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
                exp_id=`parse_property.py -f "'${bam_set[$ix]}'" --project "${DX_PROJECT_CONTEXT_ID}" --exp_id -q`
            fi
            rep_tech=`parse_property.py -f "'${bam_set[$ix]}'" --project "${DX_PROJECT_CONTEXT_ID}" --rep_tech -q`
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
        #samtools sort -@ $nthreads -m 6G -f sofar.bam sorted.bam
        set -x
        samtools sort -@ $nthreads -f sofar.bam sorted.bam
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
    qc_stats=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n bismark_map -f ${all_reports}`
    fi
    
    echo "* Collect bam stats..."
    set -x
    samtools flagstat ${target_root}.bam > ${target_root}_flagstat.txt
    samtools stats ${target_root}.bam > ${target_root}_samstats.txt
    head -3 ${target_root}_samstats.txt
    grep ^SN ${target_root}_samstats.txt | cut -f 2- > ${target_root}_samstats_summary.txt
    set +x

    echo "* Prepare metadata..."
    reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n samtools_flagstats -f ${target_root}_flagstat.txt`
        qc_stats=`echo $qc_stats, $meta`
        reads=`qc_metrics.py -n samtools_flagstats -f ${target_root}_flagstat.txt -k total`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${target_root}_samstats_summary.txt`
        read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${target_root}_samstats_summary.txt -k "average length"`
        qc_stats=`echo $qc_stats, $meta`
    fi
    # All qc to one file per target file:
    echo "===== samtools flagstat =====" > ${target_root}_qc.txt
    cat ${target_root}_flagstat.txt     >> ${target_root}_qc.txt
    echo " "                            >> ${target_root}_qc.txt
    echo "===== samtools stats ====="   >> ${target_root}_qc.txt
    cat ${target_root}_samstats.txt     >> ${target_root}_qc.txt


    #echo "* Download and uncompress files..."
    #bam_root=`dx describe "$bismark_bam" --name | cut -d'.' -f1`
    #bam_root=${bam_root%_bismark_biorep}
    #bam_root=${bam_root%_bismark}
    #target_root="${bam_root}_bismark"
    #dx download "$bismark_bam" -o - | samtools view - > ${target_root}.sam
    #set -x

    # NOTE: Better to use sam and let extractor use more threads  
    #echo "* Decompressing biorep bam..."
    #set -x
    #samtools view ${target_root}.bam > ${target_root}.sam
    #set +x

    # NOTE: Better to use sam and let extractor use more threads, but STORAGE IS LIMITED  
    echo "* Index biorep bam..."
    set -x
    samtools index ${target_root}.bam
    set +x

    echo "* Download and uncompress index..."
    set -x
    dx download "$dme_ix" -o - | tar zxvf -
    set +x

    echo "* Analyse methylation..."
    #gzipFlag=""
    #if [ "$gzip" == "true" ]; then
    #    echo '* Adding gzip flag'
    #    gzipFlag="--gzip"
    #else
    #    echo '* No gzip flag'
    #fi

    #bismark_methylation_extractor "$gzipFlag" -s --comprehensive --cytosine_report --CX_context --ample_mem \
    #  --output /home/dnanexus/output/ --zero_based --genome_folder input ${target_root}.sam
    # TODO missing: --no_overlap/--include_overlap --ignore_XXX  --multicore $nthreads
    # NOTE: reading a bam and outputting .gz will triple the number of cores used on multi-core.
    set -x
    mkdir -p /home/dnanexus/output/
    bismark_methylation_extractor --multicore $nthreads --paired-end -s --comprehensive --cytosine_report \
        --CX_context --ample_mem --output /home/dnanexus/output/ --zero_based --genome_folder input ${target_root}.bam
    rm ${target_root}.bam* # STORAGE IS LIMITED 
    set +x

    echo "* Create beds..."
    # NOTE: must end in ".CX_report"
    ls -l /home/dnanexus/output/
    set -x
    cxrepo-bed.py -o /home/dnanexus/output /home/dnanexus/output/${target_root}.CX_report.txt
    set +x
    ls -l /home/dnanexus/output/
    set -x
    mv /home/dnanexus/output/CG_${target_root}.CX_report.txt  ${target_root}_CpG.bed
    mv /home/dnanexus/output/CHG_${target_root}.CX_report.txt ${target_root}_CHG.bed
    mv /home/dnanexus/output/CHH_${target_root}.CX_report.txt ${target_root}_CHH.bed
    mv /home/dnanexus/output/${target_root}.M-bias.txt ${target_root}_mbias_report.txt
    set +x
    ls -l 

    echo "* Convert to BigBed..."
    set -x
    bedToBigBed ${target_root}_CpG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 input/chrom.sizes ${target_root}_CpG.bb
    bedToBigBed ${target_root}_CHG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 input/chrom.sizes ${target_root}_CHG.bb
    bedToBigBed ${target_root}_CHH.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 input/chrom.sizes ${target_root}_CHH.bb
    gzip *.bed
    set +x
    ls -l 

    echo "* Prepare metadata..."
    ## TODO: Figure out what the metadata is
    qc_stats=''
    #if [ -f /usr/bin/qc_metrics.py ]; then
    #    qc_stats=`qc_metrics.py -n bismark_map -f ${target_root}_map_report.txt`
    #fi
    
    echo "* Uploading files..."
    #bam_biorep=$(dx upload ${target_root}.bam --details "{ $qc_stats }" --property SW="$versions" \
    #                                           --property reads="$reads" --property read_length="$read_len" --brief)
    bam_biorep_qc=$(dx upload ${target_root}_qc.txt --details "{ $qc_stats }" --property SW="$versions" --brief)
    map_biorep=$(dx upload ${target_root}_map_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)

    #dx-jobutil-add-output bam_biorep "$bam_biorep" --class=file
    dx-jobutil-add-output bam_biorep_qc "$bam_biorep_qc" --class=file
    dx-jobutil-add-output map_biorep "$map_biorep" --class=file

    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string
    
    CpG_bed=$(dx upload ${target_root}_CpG.bed.gz --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHG_bed=$(dx upload ${target_root}_CHG.bed.gz --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHH_bed=$(dx upload ${target_root}_CHH.bed.gz --details "{ $qc_stats }" --property SW="$versions" --brief)
    CpG_bb=$(dx upload ${target_root}_CpG.bb --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHG_bb=$(dx upload ${target_root}_CHG.bb --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHH_bb=$(dx upload ${target_root}_CHH.bb --details "{ $qc_stats }" --property SW="$versions" --brief)

    mbias_report=$(dx upload ${target_root}_mbias_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output CpG_bed "$CpG_bed" --class=file
    dx-jobutil-add-output CHG_bed "$CHG_bed" --class=file
    dx-jobutil-add-output CHH_bed "$CHH_bed" --class=file
    dx-jobutil-add-output CpG_bb "$CpG_bb" --class=file
    dx-jobutil-add-output CHG_bb "$CHG_bb" --class=file
    dx-jobutil-add-output CHH_bb "$CHH_bb" --class=file
    dx-jobutil-add-output mbias_report "$mbias_report" --class=file

    #dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
