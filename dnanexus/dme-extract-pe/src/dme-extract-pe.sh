#!/bin/bash
# dme-extract-pe.sh - WGBS ENCODE Pipeline step: Extract paired-end methylation and report Whole Genome Bisulphite Analysis.

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bismark_bam: '$bismark_bam'"
    echo "* Value of dme_ix:      '$dme_ix'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of gzip:        '$gzip'"
    echo "* Value of nthreads:    '$nthreads'"

    echo "* Download and uncompress files..."
    bam_root=`dx describe "$bismark_bam" --name | cut -d'.' -f1`
    bam_root=${bam_root%_bismark_biorep}
    bam_root=${bam_root%_bismark}
    target_root="${bam_root}_bismark"
    dx download "$bismark_bam" -o - | samtools view -u - ${target_root}.sam

    dx download "$dme_ix" -o - | tar zxvf -
    dx download "$chrom_sizes" -o chrom.sizes

    echo "* Analyse methylation..."
    gzipFlag=""
    if [ "$gzip" == "true" ]; then
        echo '* Adding gzip flag'
        gzipFlag="--gzip"
    else
        echo '* No gzip flag'
    fi

    #bismark_methylation_extractor "$gzipFlag" -s --comprehensive --cytosine_report --CX_context --ample_mem \
    #  --output /home/dnanexus/output/ --zero_based --genome_folder input ${target_root}.sam
    # TODO missing: --no_overlap/--include_overlap --ignore_XXX  --multicore $nthreads
    set -x
    bismark_methylation_extractor "$gzipFlag" --paired-end -s --comprehensive --cytosine_report --CX_context --ample_mem \
      --output /home/dnanexus/output/ --zero_based --genome_folder input ${target_root}.sam
    set +x

    echo "* Create beds..."
    # NOTE: must end in ".CX_report"
    set -x
    cxrepo-bed.py -o /home/dnanexus/output ${target_root}.CX_report.txt
    set +x
    echo `ls /home/dnanexus/output`
    set -x
    mv /home/dnanexus/output/CG_${target_root}.CX_report.txt ${target_root}_CG.bed
    mv /home/dnanexus/output/CHG_${target_root}.CX_report.txt ${target_root}_CHG.bed
    mv /home/dnanexus/output/CHH_${target_root}.CX_report.txt ${target_root}_CHH.bed
    set +x

    echo "* Convert to BigBed..."
    set -x
    bedToBigBed ${target_root}_CG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes ${target_root}_CG.bb
    bedToBigBed ${target_root}_CHG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes ${target_root}_CHG.bb
    bedToBigBed ${target_root}_CHH.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes ${target_root}_CHH.bb
    set +x

    echo "* Prepare metadata..."
    # TODO: Warning, map_report may be output of mapping step and *E_ is probably SE_ and PE_
    set -x
    gzip *.bed
    cat output/*E_report.txt > ${target_root}_map_report.txt
    set +x
    # TODO: Do we need another map report here?  Was this done in merge?
    qc_stats=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n bismark_map -f ${target_root}_map_report.txt`
    fi
    set -x
    mv /home/dnanexus/output/${target_root}.M-bias.txt ${target_root}_mbias_report.txt
    gzip *.bed
    set +x
    
    
    echo "* Uploading files..."
    CG_bed=$(dx upload ${target_root}_CG.bed.gz --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHG_bed=$(dx upload ${target_root}_CHG.bed.gz --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHH_bed=$(dx upload ${target_root}_CHH.bedgz --details "{ $qc_stats }" --property SW="$versions" --brief)

    CG_bb=$(dx upload ${target_root}_CG.bb --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHG_bb=$(dx upload ${target_root}_CHG.bb --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHH_bb=$(dx upload ${target_root}_CHH.bb --details "{ $qc_stats }" --property SW="$versions" --brief)

    map_report=$(dx upload ${target_root}_map_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)
    M_bias_report=$(dx upload ${target_root}_mbias_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output CG_bed "$CG_bed" --class=file
    dx-jobutil-add-output CHG_bed "$CHG_bed" --class=file
    dx-jobutil-add-output CHH_bed "$CHH_bed" --class=file
    dx-jobutil-add-output CG_bb "$CG_bb" --class=file
    dx-jobutil-add-output CHG_bb "$CHG_bb" --class=file
    dx-jobutil-add-output CHH_bb "$CHH_bb" --class=file
    dx-jobutil-add-output map_report "$map_report" --class=file
    dx-jobutil-add-output mbias_report "$mbias_report" --class=file
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
