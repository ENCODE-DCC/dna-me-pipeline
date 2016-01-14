#!/bin/bash
# dme-bg-to-signal.sh - WGBS ENCODE Pipeline step: Convert methylation bedGraph to bigWig signal file.

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bg_gz:       '$bg_gz'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download and uncompress files..."
    bg_fn=`dx describe "$bg_gz" --name`
    target_root=${bg_fn%.bedGraph.gz}
    echo "* BedGraph file: '"$target_root".bg.gz'"
    dx download "$bg_gz" -o ${target_root}.bg.gz
    set -x
    gunzip ${target_root}.bg.gz
    set +x

    dx download "$chrom_sizes" -o chrom.sizes

    echo "* Parse metadata..."
    qc_stats=''
    if [ -f /usr/bin/parse_property.py ]; then
        qc_stats=`parse_property.py -f "'${bg_gz}'" --details`
    fi

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dname_bg_to_signal.sh ${target_root}.bg chrom.sizes
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    
    echo "* Uploading files..."
    signal=$(dx upload ${target_root}.bw --details "$qc_stats" --property SW="$versions" --brief)
    dx-jobutil-add-output signal "$signal" --class=file
    
    dx-jobutil-add-output metadata "$qc_stats" --class=string

    echo "* Finished."
}
