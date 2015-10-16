#!/bin/bash
# dme-extract-se.sh - WGBS ENCODE Pipeline step: Extract single-ended methylation and report Whole Genome Bisulphite Analysis.

main() {
    # ====================== Experimental...
    echo "* Installing macs2 and dependencies (numpy)..." 2>&1 | tee -a install.log
    #wget http://sourceforge.net/projects/numpy/files/NumPy/1.6.2/numpy-1.6.2.tar.gz/download -O numpy.tgz 2>&1 | tee -a install.log
    #mkdir numpy
    #tar -xzf numpy.tgz -C numpy --strip-components=1
    #cd numpy
    #python2.7 setup.py build --fcompiler=gnu95 2>&1 | tee -a install.log
    #sudo python2.7 setup.py install 2>&1 | tee -a install.log
    #pip install scipy
    #pip install numpy
    #pip install cutadapt
    #pip install --user --upgrade cutadapt
    #  "execDepends": [
	#    {"name": "python-numpy"},
	#    {"name": "scipy", "package_manager": "pip"},
	#    {"name": "cutadapt", "package_manager": "pip"}
    #  ],
    #git clone https://yphe@bitbucket.org/schultzmattd/methylpy.git
    #export PYTHONPATH=./methylpy/:$PYTHONPATH
    # ... Experimental =====================

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bismark_bam: '$bismark_bam'"
    echo "* Value of dme_ix:      '$dme_ix'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of nthreads:    '$nthreads'"

    echo "* Download and uncompress files..."
    bam_root=`dx describe "$bismark_bam" --name | cut -d'.' -f1`
    bam_root=${bam_root%_bismark_biorep}
    bam_root=${bam_root%_bismark}
    target_root="${bam_root}_bismark"
    set -x
    #mkdir /home/dnanexus/output/
    dx download "$bismark_bam" -o - | samtools view - > ${target_root}.sam

    dx download "$dme_ix" -o - | tar -zxvf -
    dx download "$chrom_sizes" -o chrom.sizes
    set +x

    ls -l
    
    echo "* Analyse methylation..."
    #gzipFlag=""
    #if [ "$gzip" == "true" ]; then
    #    echo '* Adding gzip flag'
    #    gzipFlag="--gzip"
    #else
    #    echo '* No gzip flag'
    #fi

    #bismark_methylation_extractor "$gzipFlag" -s --comprehensive --cytosine_report --CX_context --ample_mem \
    #  --output /home/dnanexus/output/ --zero_based --genome_folder input ${target_root}.sam
    # TODO: missing: --no_overlap/--include_overlap --ignore_XXX  --bedGraph
    # NOTE: reading a bam and outputting .gz will triple the number of cores used on multi-core.
    set -x
    mkdir -p /home/dnanexus/output/
    bismark_methylation_extractor --multicore $nthreads --single-end -s --comprehensive --cytosine_report \
        --CX_context --ample_mem --output /home/dnanexus/output/ --zero_based --genome_folder input ${target_root}.sam
    set +x

    echo "* Create beds..."
    # NOTE: must end in ".CX_report"
    ls -l /home/dnanexus/output/
    set -x
    cxrepo-bed.py -o /home/dnanexus/output /home/dnanexus/output/${target_root}.CX_report.txt
    set +x
    ls -l /home/dnanexus/output/
    set -x
    mv /home/dnanexus/output/CG_${target_root}.CX_report.txt ${target_root}_CG.bed
    mv /home/dnanexus/output/CHG_${target_root}.CX_report.txt ${target_root}_CHG.bed
    mv /home/dnanexus/output/CHH_${target_root}.CX_report.txt ${target_root}_CHH.bed
    mv /home/dnanexus/output/${target_root}.M-bias.txt ${target_root}_mbias_report.txt
    set +x
    ls -l 

    echo "* Convert to BigBed..."
    set -x
    bedToBigBed ${target_root}_CG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes ${target_root}_CG.bb
    bedToBigBed ${target_root}_CHG.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes ${target_root}_CHG.bb
    bedToBigBed ${target_root}_CHH.bed -as=/opt/data/as/bedMethyl.as -type=bed9+2 chrom.sizes ${target_root}_CHH.bb
    gzip *.bed
    set +x
    ls -l 

    echo "* Prepare metadata..."
    ## TODO: Figure out what the metadata is
    qc_stats=''
    #if [ -f /usr/bin/qc_metrics.py ]; then
    #    qc_stats=`qc_metrics.py -n bismark_map -f ${target_root}_map_report.txt`
    #fi    
    
    echo "* Uploading files..."
    CG_bed=$(dx upload ${target_root}_CG.bed.gz --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHG_bed=$(dx upload ${target_root}_CHG.bed.gz --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHH_bed=$(dx upload ${target_root}_CHH.bed.gz --details "{ $qc_stats }" --property SW="$versions" --brief)

    CG_bb=$(dx upload ${target_root}_CG.bb --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHG_bb=$(dx upload ${target_root}_CHG.bb --details "{ $qc_stats }" --property SW="$versions" --brief)
    CHH_bb=$(dx upload ${target_root}_CHH.bb --details "{ $qc_stats }" --property SW="$versions" --brief)

    M_bias_report=$(dx upload ${target_root}_mbias_report.txt --details "{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output CG_bed "$CG_bed" --class=file
    dx-jobutil-add-output CHG_bed "$CHG_bed" --class=file
    dx-jobutil-add-output CHH_bed "$CHH_bed" --class=file
    dx-jobutil-add-output CG_bb "$CG_bb" --class=file
    dx-jobutil-add-output CHG_bb "$CHG_bb" --class=file
    dx-jobutil-add-output CHH_bb "$CHH_bb" --class=file
    dx-jobutil-add-output mbias_report "$mbias_report" --class=file
    dx-jobutil-add-output metadata "{ $versions }" --class=string

    echo "* Finished."
}
