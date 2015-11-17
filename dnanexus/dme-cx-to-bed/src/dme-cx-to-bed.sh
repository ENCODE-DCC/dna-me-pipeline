#!/bin/bash
# dme-cx-to-bed.sh - WGBS ENCODE Pipeline step: Convert methylation context report to bed and bigBed files.

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of cx_report:   '$cx_report'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    cx_fn=`dx describe "$cx_report" --name`
    target_root=${cx_fn%.CX_report.txt.gz}
    echo "* Context file: '"$target_root".CX_report.txt.gz'"
    mkdir -p output
    dx download "$cx_report" -o output/${target_root}.CX_report.txt.gz

    dx download "$chrom_sizes" -o chrom.sizes
    
    echo "* Parse metadata..."
    qc_stats=''
    reads=0
    #read_len=0
    if [ -f /usr/bin/parse_property.py ]; then
        qc_stats=`parse_property.py -f "'${cx_report}'" --details`
    fi

    echo "* Create beds..."
    set -x
    gunzip output/${target_root}.CX_report.txt.gz
    cxrepo-bed.py -o output/ output/${target_root}.CX_report.txt
    set +x
    ls -l output/
    set -x
    mv output/CG_${target_root}.CX_report.txt  ${target_root}_CpG.bed
    mv output/CHG_${target_root}.CX_report.txt ${target_root}_CHG.bed
    mv output/CHH_${target_root}.CX_report.txt ${target_root}_CHH.bed
    set +x
    ls -l 

    echo "* Convert to BigBed..."
    set -x
    bedToBigBed ${target_root}_CHH.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes ${target_root}_CHH.bb
    pigz ${target_root}_CHH.bed
    bedToBigBed ${target_root}_CHG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes ${target_root}_CHG.bb
    pigz ${target_root}_CHG.bed
    bedToBigBed ${target_root}_CpG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes ${target_root}_CpG.bb
    pigz ${target_root}_CpG.bed
    set +x
    echo "* Check storage..."
    ls -l 
    df -k .

    echo "* Uploading files..."
    CpG_bed=$(dx upload ${target_root}_CpG.bed.gz --details "$qc_stats" --property SW="$versions" --brief)
    CHG_bed=$(dx upload ${target_root}_CHG.bed.gz --details "$qc_stats" --property SW="$versions" --brief)
    CHH_bed=$(dx upload ${target_root}_CHH.bed.gz --details "$qc_stats" --property SW="$versions" --brief)
    CpG_bb=$(dx upload ${target_root}_CpG.bb --details "$qc_stats" --property SW="$versions" --brief)
    CHG_bb=$(dx upload ${target_root}_CHG.bb --details "$qc_stats" --property SW="$versions" --brief)
    CHH_bb=$(dx upload ${target_root}_CHH.bb --details "$qc_stats" --property SW="$versions" --brief)

    dx-jobutil-add-output CpG_bed "$CpG_bed" --class=file
    dx-jobutil-add-output CHG_bed "$CHG_bed" --class=file
    dx-jobutil-add-output CHH_bed "$CHH_bed" --class=file
    dx-jobutil-add-output CpG_bb "$CpG_bb" --class=file
    dx-jobutil-add-output CHG_bb "$CHG_bb" --class=file
    dx-jobutil-add-output CHH_bb "$CHH_bb" --class=file

    dx-jobutil-add-output metadata "$qc_stats" --class=string

    echo "* Finished."
}
