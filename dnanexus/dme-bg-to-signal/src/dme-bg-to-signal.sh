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

    echo "* Download files..."
    bg_fn=`dx describe "$bg_gz" --name`
    target_root=${bg_fn%.bedGraph.gz}
    echo "* BedGraph file: '"$target_root".bedGraph.gz'"
    mkdir -p output
    dx download "$bg_gz" -o output/${target_root}.bedGraph.gz

    dx download "$chrom_sizes" -o chrom.sizes

    echo "* Parse metadata..."
    qc_stats=''
    if [ -f /usr/bin/parse_property.py ]; then
        qc_stats=`parse_property.py -f "'${bg_gz}'" --details`
    fi

    echo "* Convert to signal bedGraph to bigWig..."
    # NOTE: Not sure if we want signal
    ls -l output/
    set -x
    gunzip output/${target_root}.bedGraph.gz
    bedGraphToBigWig output/${target_root}.bedGraph chrom.sizes ${target_root}.bw
    rm -f output/${target_root}.bedGraph
    set +x
    echo "* Check storage..."
    ls -l 
    df -k .

    echo "* Uploading files..."
    signal=$(dx upload ${target_root}.bw --details "$qc_stats" --property SW="$versions" --brief)
    dx-jobutil-add-output signal "$signal" --class=file
    
    dx-jobutil-add-output metadata "$qc_stats" --class=string

    echo "* Finished."
}
