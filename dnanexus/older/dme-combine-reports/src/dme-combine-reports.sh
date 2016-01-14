#!/bin/bash
# dme-align-se.sh - align se reads with bismark/bowtie1

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_set:        '$map_report_genome'"
    echo "* Value of map_report_set: '$map_report_lambda'"

    echo "* Downloading files..."
    file_root=`dx describe "$map_report_genome" --name`
    file_root=${file_root%.fq_bismark_map_report.txt}
    dx download "$map_report_genome" -o genome_bismark_map_report.txt
    dx download "$map_report_lambda" -o lambda_bismark_map_report.txt

    combined_file=${file_root}_techrep_bismark_map_report.txt
    echo "* Merged map_report file will be: '${combined_file}'"

    echo "* Prepare metadata..."
    echo "===== bismark reference =====" > ${combined_file}
    cat genome_bismark_map_report.txt   >> ${combined_file}
    echo " "                            >> ${combined_file}
    echo "===== bismark lambda ====="   >> ${combined_file}
    cat lambda_bismark_map_report.txt   >> ${combined_file}

    qc_stats=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n bismark_ref -f ${combined_file}`
    fi

    echo "* Upload results..."
    map_techrep=$(dx upload ${combined_file} --details "{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output map_techrep "$map_techrep" --class=file
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
 }
