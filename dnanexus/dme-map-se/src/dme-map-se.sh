#!/bin/bash
# dme-map-se.sh

set -x
set +e

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of trimmed_reads: '$trimmed_reads'"
    echo "* Value of dme_ix:        '$dme_ix'"
    #echo "* Value of min_insert:     $min_insert"
    #echo "* Value of max_insert:     $max_insert"

    echo "* Download files..."
    reads_root=`dx describe "$trimmed_reads" --name | cut -d'.' -f1`
    dx download "$trimmed_reads" -o - | gunzip > "$reads_root".fq
    dx download "$dme_ix" -o - | tar zxvf -

    bam_root="${reads_root}_bismark"
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
    bismark -n 1 -l 28 -output_dir output --temp_dir output input "$reads_root".fq
    set +x

    echo "* Compressing with samtools..."
    set -x
    samtools -Sb -@ ${nthreads} output/${reads_root}.fq_bismark_se.sam ${bam_root}.bam
    set +x


    echo "* Upload results..."
    ls -l /home/dnanexus/output
    #echo `ls /home/dnanexus/output`
    #outfile="$read_fn".mapped_methylseq.tgz
    #tar zcvf $outfile output
    #mapped_files=$(dx upload /home/dnanexus/$outfile --brief)
    bismark_bam=$(dx upload /home/dnanexus/${bam_root}.bam --property SW="$versions" --brief)

    dx-jobutil-add-output bismark_bam "$bismark_bam" --class=file

    echo "* Finished."
 }
