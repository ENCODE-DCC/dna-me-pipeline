#!/bin/bash
# dme-extract-se.sh - WGBS ENCODE Pipeline step: Extract single-ended methylation and report Whole Genome Bisulphite Analysis.

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
        set -x
        samtools sort -@ $nthreads -m 6G -f sofar.bam sorted.bam
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
    
    echo "* Collect bam stats..."
    set -x
    samtools flagstat ${target_root}.bam > ${target_root}_flagstat.txt
    set +x
    # NOTE: samtools stats may take longer than it is worth 
    #samtools stats ${target_root}.bam > ${target_root}_samstats.txt
    #head -3 ${target_root}_samstats.txt
    #grep ^SN ${target_root}_samstats.txt | cut -f 2- > ${target_root}_samstats_summary.txt

    echo "* Prepare metadata..."
    qc_stats=''
    reads=0
    #read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n bismark_map -f ${all_reports}`
        meta=`qc_metrics.py -n samtools_flagstats -f ${target_root}_flagstat.txt`
        qc_stats=`echo $qc_stats, $meta`
        reads=`qc_metrics.py -n samtools_flagstats -f ${target_root}_flagstat.txt -k total`
        #meta=`qc_metrics.py -n samtools_stats -d ':' -f ${target_root}_samstats_summary.txt`
        #read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${target_root}_samstats_summary.txt -k "average length"`
        #qc_stats=`echo $qc_stats, $meta`
    fi
    # All qc to one file per target file:
    echo "===== samtools flagstat =====" > ${target_root}_qc.txt
    cat ${target_root}_flagstat.txt     >> ${target_root}_qc.txt
    #echo " "                            >> ${target_root}_qc.txt
    #echo "===== samtools stats ====="   >> ${target_root}_qc.txt
    #cat ${target_root}_samstats.txt     >> ${target_root}_qc.txt

    # NOTE: Better to use sam and let extractor use more threads, but this takes up precious storage
    alignments_file=${target_root}.bam
    ncores=$nthreads
    # TODO: What bam_size constitutes too large for sam?  93,953,130,496 is fine!
    bam_size=`ls -go ${target_root}.bam | awk '{print $3}'`
    if [ "$uncompress_bam" == "true" ] && [ $bam_size -lt 400000000000 ]; then
        alignments_file=${target_root}.sam
        echo "* Decompressing biorep bam (size: $bam_size)..."
        set -x
        samtools view ${target_root}.bam > ${alignments_file}
        rm ${target_root}.bam
        set +x
    else
        ncores=`expr $nthreads / 2`
        if [ "$uncompress_bam" != "true" ]; then
            echo "* Using compressed biorep bam (size: $bam_size) and $ncores cores..."
        else
            echo "* Using compressed bam and $ncores cores because bam_size: $bam_size exceeds limit."
        fi
    fi

    echo "* Download and uncompress index..."
    set -x
    dx download "$dme_ix" -o - | tar -zxf -
    set +x

    echo "* Analyse methylation in ${alignments_file} and using $ncores threads..."
    # TODO: missing: --no_overlap/--include_overlap --ignore_XXX
    # NOTE: reading a bam and outputting .gz will triple the number of cores used on multi-core.
    set -x
    mkdir -p /home/dnanexus/output/
    bismark_methylation_extractor --multicore $ncores --single-end -s --comprehensive --cytosine_report \
        --CX_context --ample_mem --output output/ --zero_based --genome_folder input ${alignments_file}
    rm -f ${alignments_file} # STORAGE IS LIMITED
    set +x

    echo "* More metadata..."
    ls -l output/*_splitting_report.txt
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n bismark_extract -f output/*_splitting_report.txt`
        qc_stats=`echo $qc_stats, $meta`
    fi
    echo " "                                                             >> ${target_root}_qc.txt
    echo "===== bismark_methylation_extractor: splitting_report ====="   >> ${target_root}_qc.txt
    cat output/${target_root}_splitting_report.txt                       >> ${target_root}_qc.txt

    echo "* Convert to signal bedGraph to bigWig..."
    # NOTE: Not sure if we want signal
    ls -l output/
    set -x
    rm -f output/*_context_${target_root}.txt # STORAGE IS LIMITED
    rm -f output/${target_root}.bam_splitting_report.txt
    rm -f output/${target_root}.bedGraph.gz.bismark.zero.cov
    rm -f output/${target_root}.bismark.cov.gz
    gunzip output/${target_root}.bedGraph.gz
    bedGraphToBigWig output/${target_root}.bedGraph input/chrom.sizes ${target_root}.bw
    rm -f output/${target_root}.bedGraph
    set +x
    
    echo "* Create beds..."
    set -x
    cxrepo-bed.py -o output/ output/${target_root}.CX_report.txt
    set +x
    ls -l output/
    set -x
    mv output/CG_${target_root}.CX_report.txt  ${target_root}_CpG.bed
    mv output/CHG_${target_root}.CX_report.txt ${target_root}_CHG.bed
    mv output/CHH_${target_root}.CX_report.txt ${target_root}_CHH.bed
    mv output/${target_root}.M-bias.txt ${target_root}_mbias_report.txt
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
    echo "* Check storage..."
    df -k .

    echo "* Uploading files..."
    # NOTE: Not saving merged bam
    #bam_biorep=$(dx upload ${target_root}.bam --details "{ $qc_stats }" --property SW="$versions" \
    #                                           --property reads="$reads" --property read_length="$read_len" --brief)
    bam_biorep_qc=$(dx upload ${target_root}_qc.txt --details "{ $qc_stats }" --property SW="$versions" --brief)
    map_biorep=$(dx upload ${target_root}_map_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)

    #dx-jobutil-add-output bam_biorep "$bam_biorep" --class=file
    dx-jobutil-add-output bam_biorep_qc "$bam_biorep_qc" --class=file
    dx-jobutil-add-output map_biorep "$map_biorep" --class=file

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

    # NOTE: Not sure if we want signal
    signal=$(dx upload ${target_root}.bw --details "{ $qc_stats }" --property SW="$versions" --brief)
    dx-jobutil-add-output signal "$signal" --class=file
    
    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
