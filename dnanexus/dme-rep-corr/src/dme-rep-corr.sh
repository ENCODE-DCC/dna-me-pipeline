#!/bin/bash
# dme-rep-corr.sh - WGBS ENCODE Pipeline step: Correlates replicate CpG methylation files.

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bg_gz:       '$CpG_A'"
    echo "* Value of bg_gz:       '$CpG_B'"

    echo "* Download files..."
    A_file=`dx describe "$CpG_A" --name`
    A_root=${A_file%.gz}
    A_root=${A_root%.bed}
    A_root=${A_root%_bismark_biorep_CpG}
    A_root=${A_root%_bismark_CpG}
    echo "* A file: '$A_file'"
    dx download "$CpG_A" -o ${A_root}_A.bed.gz

    echo "* Download files..."
    B_file=`dx describe "$CpG_B" --name`
    B_root=${B_file%.gz}
    B_root=${B_root%.bed}
    B_root=${B_root%_bismark_biorep_CpG}
    B_root=${B_root%_bismark_CpG}
    echo "* B file: '$B_file'"
    dx download "$CpG_B" -o ${B_root}_B.bed.gz
    
    target_root="${A_root}_${B_root}"
    if [ -f /usr/bin/parse_property.py ]; then
        exp_id=`parse_property.py -f "'$CpG_A'" --project "${DX_PROJECT_CONTEXT_ID}" --exp_id --quiet`
        rep_A=`parse_property.py -f "'$CpG_A'" --project "${DX_PROJECT_CONTEXT_ID}" --rep_tech --quiet`
        rep_B=`parse_property.py -f "'$CpG_B'" --project "${DX_PROJECT_CONTEXT_ID}" --rep_tech --quiet`
        if [ "$exp_id" != "" ] && [ "$rep_A" != "" ] && [ "$rep_B" != "" ]; then
            target_root="${exp_id}_${rep_A}_${rep_B}"
        fi
    fi

    echo "* Parse metadata..."
    qc_stats=''
    if [ -f /usr/bin/parse_property.py ]; then
        qc_stats=`parse_property.py -f "'${bg_gz}'" --details`
    fi

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dname_bed_corr.sh ${A_root}_A.bed.gz ${B_root}_B.bed.gz "CpG" $target_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    target_root="${target_root}_CpG_corr"
    
    echo "* Prepare metadata..."
    qc_stats=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n bedmethyl_corr -t vertical -d = -f ${target_root}.txt`
    fi

    echo "* Uploading files..."
    CpG_corr=$(dx upload ${target_root}.txt --details "{ $qc_stats }" --property SW="$versions" --brief)
    dx-jobutil-add-output CpG_corr "$CpG_corr" --class=file
    
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
