#!/bin/bash

### Instruction: This script is being used by other script in bismark_pipeline dir call_mergeUnsorted_dedup_files_for_methExtraction.sh, and this script would merge the unsorted bam files, 
### and deduplicate the unsorted_merged bam files, which would be further process for the methylation calling. More importantly, since methylation extraction doesn't want 
### any positional sorting to be done, so you would want to merge all the unsorted bam files, output of bismark alignment, in to one merged bam file for methylation call later on.
### Basically, you would be skipping the sorting and indexing step here, which could save considerable computation time too.	

export RUN_PATH=`pwd`

for LIB in ${LIST}; do
		
	if [[ ! -d $OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file ]]; then
		mkdir $OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file
	fi

	BAM_COUNT=$( ls $BAM_DIR/$LIB/bam_files/*.bam | wc -l)
	UNSORTED_BAM_DIR=$BAM_DIR/$LIB/unsortedButMerged_ForBismark_file
	
	### Merge the unsorted bam files:
	if [[ $BAM_COUNT > 1 ]]; then
                BAM_LIST=$( ls $BAM_DIR/$LIB/bam_files/*.bam | tr '\n' ' ' )
                $SAMTOOLS_PATH/samtools merge -nf $OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file/${LIB}_unsorted_merged.bam $BAM_LIST

        else if [[ $BAM_COUNT = 1 ]]; then
                cp $BAM_DIR/$LIB/bam_files/*.bam $OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file/${LIB}_unsorted_merged.bam
                
		fi      
        fi     
	
	### Run the deduplication, and remove the pcr duplicates from unsorted_bam_files (i.e the sequence aligning to the same genomic positions).
	$BISMARK_PATH/deduplicate_bismark -p --bam $OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file/${LIB}_unsorted_merged.bam

done


for LIB in ${LIST}; do
	BAM_PATH=$INPUT_DIR/$LIB/unsortedButMerged_ForBismark_file
	
	if [[ ! -d $OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file ]]; then
		mkdir $OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file
	fi

	if [[ ! -d $OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file/methylation_extraction ]]; then	
		mkdir $OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file/methylation_extraction
	fi
	
	### Set the output dir to retain all the methylation called file:
	METH_OUTPUT_DIR=$OUTPUT_DIR/$LIB/unsortedButMerged_ForBismark_file/methylation_extraction
	echo "Processing $LIB"
	echo "Using BAM files in $BAM_PATH"
	echo "Using $METH_OUTPUT_DIR for output"
	
	### BAM files should be unsorted, and Default = CpG context only; else use --CX_context for all CpG info.
	for bam_file in $(/bin/ls $BAM_PATH/*_unsorted_merged.deduplicated.bam); do
		$BISMARK_PATH/bismark_methylation_extractor -p --no_overlap -o $METH_OUTPUT_DIR --comprehensive --report --bedGraph --genome_folder $GENOME_PATH $bam_file	
		echo ""
	done
done

