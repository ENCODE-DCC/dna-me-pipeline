#!/bin/bash

### Instruction: This script finds the pair for all those splitted fastqs to make it paired-end. And further runs trim-galore, and bismark on each paired-end splitted fastqs on different 
### nodes for parallel job completion using trim_galore_bismark_alignment.sh script. Moreover, the parameter for job run has been designed in such a way that the job submission of this 
### trim_galore_bismark_alignment.sh to the morgan cluster is dependent on completion of all fastq splitting executed by prior fastq splitting based scripts.
  

### LIB_LIST environment variable, containing the list of SL# library with the wgbs fastqs, created in prior script imported:
for LIB in ${LIB_LIST}; do
	DIR_PATH=$OUTPUT_DIR/$LIB
	echo "Processing ${LIB}"
	echo "Using FASTQ files in ${DIR_PATH}"
	
	### Create all sub-directories in main dir:
	if [[ ! -d $OUTPUT_DIR/${LIB} ]]; then
		mkdir $OUTPUT_DIR/${LIB}
	fi
	
	if [[ ! -d $OUTPUT_DIR/${LIB}/temp_dir ]]; then
		mkdir $OUTPUT_DIR/${LIB}/temp_dir
	fi
	
	if [[ ! -d $OUTPUT_DIR/${LIB}/bam_files ]]; then
		mkdir $OUTPUT_DIR/${LIB}/bam_files
	fi
	
	if [[ ! -d $OUTPUT_DIR/${LIB}/log_files ]]; then
		mkdir $OUTPUT_DIR/${LIB}/log_files
	fi
	
	echo "Using $OUTPUT_DIR for output, and all sub directories created.."
	export LOG_DIR=$OUTPUT_DIR/${LIB}/log_files ### Env var for later usage to store all log files
	BAM_OUTPUT_DIR=$OUTPUT_DIR/${LIB}/bam_files

	### Find the pair for each splittted fastq files to make it paired-end and feed it into the bismark alignment pipeline for further processing using trim_galore_bismark_alignment.sh script:
	### The fastq considered here would be of the format instance like : C6RLPANXX_s4_1_7bp_Index_7_SL83766.00 
	for fastq in $( ls $DIR_PATH/*); do
			#FASTQ_ID=$(basename $fastq | sed 's/\./_/g'| awk -F"_" '{print $1"_"$2"_"$NF}' )
			SEQUENCER_ID=$(basename $fastq | cut -f 1 -d '_')
			FLOWLANE_ID=$(basename $fastq | cut -f 2 -d '_')
			READ_1=$(basename $fastq | cut -f 3 -d '_')
			SPLIT_ID=$(basename $fastq | awk -F"." '{print $NF}')


		for fastq2 in $( ls $DIR_PATH/*); do				
			#FASTQ_ID_SEARCH=$(basename $fastq2 | sed 's/\./_/g'| awk -F"_" '{print $1"_"$2"_"$NF}' )
			SEQUENCER_ID_SEARCH=$(basename $fastq2 | cut -f 1 -d '_') 
			FLOWLANE_ID_SEARCH=$(basename $fastq2 | cut -f 2 -d '_')
			READ_SEARCH=$(basename $fastq2 | cut -f 3 -d '_')
			SPLIT_ID_SEARCH=$(basename $fastq2 | awk -F"." '{print $NF}')

			if [[ $READ_1 = "1" && $SEQUENCER_ID = $SEQUENCER_ID_SEARCH && $READ_SEARCH = "2" && $SPLIT_ID = $SPLIT_ID_SEARCH && $FLOWLANE_ID = $FLOWLANE_ID_SEARCH ]]; then
				FASTQ_READ_1=$fastq
				FASTQ_READ_2=$fastq2
				LOG_BASE=$(basename $FASTQ)
				JOB_NAME="Runnning bismark"
				### Will take 8 cores of the morgan cluster for each bismark alignment instances:
				bsub -We 72:00 -q night -n $CORE_NUM -R span[hosts=1] -J $JOB_NAME -o $LOG_DIR/${LIB}_bismarkrun_log.out ./trim_galore_bismark_alignment.sh $FASTQ_READ_1 $FASTQ_READ_2 $BAM_OUTPUT_DIR 
			fi
		done
	done
done






