#!/bin/bash
# dme-map-pe.sh

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of trimmed_reads1: '$trimmed_reads1'"
    echo "* Value of trimmed_reads2: '$trimmed_reads2'"
    echo "* Value of dme_ix:         '$dme_ix'"
    echo "* Value of min_insert:      $min_insert"
    echo "* Value of max_insert:      $max_insert"

    echo "* Download files..."
    reads1_root=`dx describe "$trimmed_reads1" --name | cut -d'.' -f1`
    reads2_root=`dx describe "$trimmed_reads2" --name | cut -d'.' -f1`
    dx download "$trimmed_reads1" -o - | gunzip > "$reads1_root".fq
    dx download "$trimmed_reads2" -o - | gunzip > "$reads2_root".fq
    dx download "$dme_ix" -o - | tar zxvf -

    bam_root="${reads1_root}_${reads2_root}_bismark"
    # Try to simplify the names
    if [ -f /usr/bin/parse_property.py ]; then
        rep_root=`parse_property.py -f "'${trimmed_reads1}'" --project "${DX_PROJECT_CONTEXT_ID}" --root_name --quiet`
        if [ "$rep_root" != "" ]; then
            bam_root="${rep_root}_bismark"
        fi
    fi
    echo "* Expect to create '${bam_root}.bam'"

    echo "* Mapping with bismark/bowtie..."
    #ls -l input
    set -x
    mkdir output
    bismark -n 1 -l 28 -output_dir output --temp_dir output input \
            -I $min_insert -X $max_insert \
            -1 ${reads1_root}.fq -2 ${reads2_root}.fq
    set +x

    echo "* Compressing with samtools..."
    set -x
    samtools -Sb -@ ${nthreads} output/${reads1_root}.fq_bismark_pe.sam ${bam_root}.bam
    set +x

    echo "* Upload results..."
    ls -l /home/dnanexus/output
    #mv output/${reads1_root}.fq_bismark_pe.sam output/${bam_root}.sam
    ## rename necessary to conserve extraction step
    ## creates PE report file which is not used.
    #outfile="$read_fn1".mapped_methylseq.tgz
    #tar zcvf $outfile output
    bismark_bam=$(dx upload /home/dnanexus/${bam_root}.bam --property SW="$versions" --brief)

    dx-jobutil-add-output bismark_bam "$bismark_bam" --class=file

    echo "* Finished."
 }
