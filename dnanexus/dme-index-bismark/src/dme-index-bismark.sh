#!/bin/bash
# dme-index-bismark.sh - index genome with bismark/bowtie1

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi
    
    echo "* Value of reference: '$reference'"
    echo "* Value of lambda:    '$lambda'"
    
    # Prefer to discover genome and gender
    source_msg="Value of"
    if [ -f /usr/bin/parse_property.py ]; then 
        genome_prop=`parse_property.py -f "$reference" -p "genome" --quiet`
        gender_prop=`parse_property.py -f "$reference" -p "gender" --quiet`
        if [ "$genome_prop" != "" ] &&  [ "$gender_prop" != "" ]; then
            genome=$genome_prop
            gender=$gender_prop
            source_msg="Discovered"
        fi
    fi
    if [ "$genome" == "" ] || [ "$gender" == "" ]; then
        echo "Reference genome and/or gender could not be determined and must be supplied as arguments."
        exit 1
    fi
    echo "* ${source_msg} genome: '$genome'"
    echo "* ${source_msg} gender: '$gender'"
    

    echo "* Download and unzip genome reference..."
    mkdir -p input/lambda
    dx download "$reference" -o - | gunzip > input/${genome}_${gender}.fa
    dx download "$lambda" -o - | gunzip > input/lambda/lambda.fa

    index_root="${genome}_${gender}_bismark_bowtie1_index"
    echo "* Expect to create '${index_root}.tgz'"
    
    echo "* Preparing/indexing ${genome}-${gender} genome..."
    set -x
    bismark_genome_preparation --bowtie1 --path_to_bowtie /usr/bin/ input | tee ref.log
    set +x
    
    echo "* Preparing/indexing lambda genome..."
    set -x
    bismark_genome_preparation --bowtie1 --path_to_bowtie /usr/bin/ input/lambda | tee lambda.log
    set +x
      
    # QC anyone?
    ref_ctot=`head -10 ref.log | grep -F "C->T" | awk '{print $2}'`
    ref_gtoa=`head -10 ref.log | grep -F "G->A" | awk '{print $2}'`
    lambda_ctot=`head -10 lambda.log | grep -F "C->T" | awk '{print $2}'`
    lambda_gtoa=`head -10 lambda.log | grep -F "G->A" | awk '{print $2}'`
    meta=`echo { \"reference\": { \"genome\": \"${genome}\", \"gender\": \"${gender}\"`
    meta=`echo ${meta}, \"C_to_T\": ${ref_ctot}, \"G_to_A\": ${ref_gtoa} }`
    meta=`echo ${meta}, \"lambda\": { \"C_to_T\": ${lambda_ctot}, \"G_to_A\": ${lambda_gtoa} } }`
    echo "* JSON metadata..."
    echo ${meta}
    echo "* ----------------"

    echo "* Archiving prepped genome..."
    ls -l input/Bisulfite_Genome/
    set -x
    tar zcvf ${index_root}.tgz input/${genome}_${gender}.fa input/Bisulfite_Genome/ input/lambda/lambda.fa input/lambda/Bisulfite_Genome/
    set +x

    echo "* Upload results..."
    dme_ix=$(dx upload ${index_root}.tgz --details "${meta}" --property genome="$genome" --property gender="$gender" \
                            --property C_to_T="$ref_ctot" --property G_to_A="$ref_gtoa" --property SW="$versions" --brief)

    dx-jobutil-add-output dme_ix "$dme_ix" --class=file
    dx-jobutil-add-output metadata "${meta}" --class=string

    echo "* Finished."
}
