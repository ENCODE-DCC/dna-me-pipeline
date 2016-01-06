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

    echo "* Download and uncompress files..."
    cx_fn=`dx describe "$cx_report" --name`
    target_root=${cx_fn%.CX_report.txt.gz}
    echo "* Context file: '"$target_root".CX_report.txt.gz'"
    dx download "$cx_report" -o ${target_root}.CX_report.txt.gz
    set -x
    gunzip ${target_root}.CX_report.txt.gz
    set +x

    dx download "$chrom_sizes" -o chrom.sizes
    
    echo "* Parse metadata..."
    qc_stats=''
    reads=0
    #read_len=0
    if [ -f /usr/bin/parse_property.py ]; then
        qc_stats=`parse_property.py -f "'${cx_report}'" --details`
    fi

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dname_cx_to_bed.sh ${target_root}.CX_report.txt chrom.sizes
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    
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
