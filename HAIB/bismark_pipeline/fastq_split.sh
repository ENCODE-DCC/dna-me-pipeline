#!/bin/bash

### Instruction: This short script is dependent, and is being called by other call_fastq_split.sh script residing in bismark_pipeline dir.


for LIB in ${LIST}; do

	### Generates the directory for each SL# inside the main output dir to retain the fastq files:	
	if [[ ! -d $OUTPUT_DIR/$LIB ]]; then
		mkdir $OUTPUT_DIR/$LIB
	fi

	### Iterates over each fastq files for each SL# to split each fastq files to 72 million lines i.e 18 million reads
	### The job of splitting each fastq from each library could also be parallelized on different nodes if the user wants, 
	### but not implemented for this instance.
	for file in $(ls $INPUT_DIR/*fastq.gz); do		
		SPLIT_NAME=$( basename $file | sed 's/fastq.gz//')
		zcat $file | split -l 72000000 -d - $OUTPUT_DIR/$LIB/$SPLIT_NAME
	done
done


