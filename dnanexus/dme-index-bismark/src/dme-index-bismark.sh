#!/bin/bash
# dme-index-bismark.sh - Creates C->T indexed genomic files with Bismark/Bowtie2 used for Whole Genome Bisulphite Analysis

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi
    
    echo "* Value of genome:     '$genome'"
    echo "* Value of lambda:     '$lambda'"
    
    # Prefer to discover assembly and gender
    source_msg="Value of"
    if [ -f /usr/bin/tool_versions.py ]; then 
        assembly_prop=`parse_property.py -f "$genome" -p "genome" --quiet`
        gender_prop=`parse_property.py -f "$genome" -p "gender" --quiet`
        if [ "$assembly_prop" != "" ] &&  [ "$gender_prop" != "" ]; then
            assembly=$assembly_prop
            gender=$gender_prop
            source_msg="Discovered"
        fi
    fi
    if [ "$assembly" == "" ] || [ "$gender" == "" ]; then
        echo "Reference assembly and/or gender could not be determined and must be supplied as arguments."
        exit 1
    fi
    echo "* ${source_msg} assembly: '$assembly'"
    echo "* ${source_msg} gender:   '$gender'"
    

    echo "* Download and unzip genome reference..."
    mkdir -p input/lambda
    dx download "$genome" -o - | gunzip > input/genome.fa
    dx download "$lambda" -o - | gunzip > input/lambda/lambda.fa

    index_root="${assembly}_${gender}_bismark_bowtie2_index"
    echo "* Expect to create '${index_root}.tgz'"
    
    echo "* Preparing/indexing ${assembly}-${gender} genome..."
    set -x
    bismark_genome_preparation --bowtie2 --path_to_bowtie /usr/bin/ input | tee ref.log
    set +x
    
    echo "* Preparing/indexing lambda genome..."
    set -x
    bismark_genome_preparation --bowtie2 --path_to_bowtie /usr/bin/ input/lambda | tee lambda.log
    set +x
    
    # QC anyone?
    ref_ctot=`grep "C\-\>T\:" ref.log | awk '{print $2}'`
    ref_gtoa=`grep "G\-\>A\:" ref.log | awk '{print $2}'`
    lambda_ctot=`grep "C\-\>T\:" lambda.log | awk '{print $2}'`
    lambda_gtoa=`grep "G\-\>A\:" lambda.log | awk '{print $2}'`
    meta="{ \"reference\": { \"assembly\": \"${assembly}\", \"gender\":\"${gender}\""
    meta="${meta}, \"C_to_T\": ${ref_ctot}, \"G_to_A\": ${ref_gtoa} }"
    meta="${meta}, \"lambda\": { \"C_to_T\": ${lambda_ctot}, \"G_to_A\": ${lambda_gtoa} } }"

    echo "* Archiving prepped genome..."
    ls -l input/Bisulfite_Genome/
    set -x
    tar zcvf ${index_root}.tgz input/genome.fa input/Bisulfite_Genome/
    set +x

    echo "* Upload results..."
    #ls -l /home/dnanexus
    dme_ix=$(dx upload ${index_root}.tgz --details "${meta}" --property genome="$assembly" --property gender="$gender" \
                            --property C_to_T="$ref_ctot" --property G_to_A="$ref_gtoa" --property SW="$versions" --brief)

    dx-jobutil-add-output dme_ix "$dme_ix" --class=file
    dx-jobutil-add-output metadata "${meta}" --class=string

    echo "* Finished."
}
