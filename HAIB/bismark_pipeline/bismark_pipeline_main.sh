#!/bin/bash

####################################################
## WGBS Pipeline for HAIB ##                       #
# Author: Surya B Chhetri                          #
#                                                  #
# Set the output dir, input dir, library list,     # 
# and paths for all the tools with this script     #   
# This script would further call on other scripts  #
# residing in the bismark_pipeline directory.      # 
####################################################


### Usage: ./call_fastq_split.sh
### Input should be provided to the variables: OUTPUT_LOC, INPUT_LOC, LIB_LIST, GENOME_PATH, BISMARK_PATH, SAMTOOLS_PATH, TRIMGALORE_PATH, BOWTIE_PATH, CORE_NUM


SOURCE_DIR="/opt/HAIB/myerslab/etc"

### Mandatory sourcing of bashrc for necessary environment variables. ###
if [ -e $SOURCE_DIR/bashrc ]; then
    . $SOURCE_DIR/bashrc
else echo "[fatal] - Could not find myerslab bashrc file. Exiting"; exit 1; fi

### Mandatory sourcing of functions to get helper functions (like call_cmd). ###
if [ -e $SOURCE_DIR/functions ]; then
    . $SOURCE_DIR/functions
else echo "[fatal] - Could not find functions file. Exiting"; exit 1; fi

### Verify we are not running on the head node. ###
if [ -z "$LSB_JOBID" ]; then log_msg fatal "Please run on a compute node. Exiting"; exit 1; fi


#####################################################
#                                                   #
#         INPUT INFORMATION REQUIRED                #
#                                                   #
#####################################################


export RUN_PATH=`pwd`

## #Provide the input dir location containing all the SL# libraries to access the fastq files of interest:
export INPUT_LOC="/gpfs/gpfs1/myerslab/data/Libraries"

### Search these SL# or library numbers containing all the wgbs fastqs under myerslab data repository:
export LIB_LIST='SL160916 SL160917 SL160918 SL160919 SL160920 SL160921'

### Set the main output dir location to retain all the splitted fastq files:
export OUTPUT_LOC="/gpfs/gpfs1/home/schhetri/wgbs_run/wgbs_split_Vth_batch"

### Set the genome path location, hg19 or grch38:
export GENOME_PATH=" /gpfs/gpfs1/home/schhetri/bismark_genome_link/bismark_genome/"

### Though bismark could be in path, set full bismark path location for consistency:
export BISMARK_PATH="/opt/bismark-0.11.1"

### Though samtools could be in path, set full samtools path location for consistency:
export SAMTOOLS_PATH="/usr/local/bin"

### Set trimgalore path location for consistency:
export TRIMGALORE_PATH="/gpfs/gpfs1/home/schhetri/Tools/trim_galore_zip"

### Set bowtie path location for consistency:
export BOWTIE_PATH="/opt/bowtie2-2.1.0"

### Set the no. of cores you want:
export CORE_NUM=8

### Create main output dir:
if [[ ! -d $OUTPUT_LOC ]]; then
	mkdir $OUTPUT_LOC
fi

export OUTPUT_DIR=$OUTPUT_LOC


#####################################################
#                                                   #
#         PIPELINE STARTS HERE                      #
#                                                   #
#####################################################


### This job calls the script call_fastq_split.sh, which has all the job submission information for fastq splitting with job name -J "Splitting of fastq":
bsub -We 24:00 -n 1 -R span[hosts=1] -J "Fastq splitting" -o $OUTPUT_DIR/fastq_split_main.out $RUN_PATH/call_fastq_split.sh

### This job won't start until the previous job with name -J "Splitting of fastq" is completed from the previous script, and upon completion eventually runs the job with name -J "Running bismark":
bsub -We 24:00 -n 1 -R span[hosts=1] -w "done('Splitting of fastq')" -J "Bismark and trim galore run for paired-end fastq" -o $OUTPUT_DIR/bismarkrun_main.out $RUN_PATH/call_trim_galore_bismark_alignment.sh  

### This job is dependent upon the completion of previous script job run with name -J "Runnning bismark":
bsub -We 24:00 -n 1 -R span[hosts=1] -w "done('Runnning bismark')" -J "Deduplication and methylation calling" -o $OUTPUT_DIR/deduplication_methylationCall_main.out $RUN_PATH/call_mergeUnsorted_dedup_files_for_methExtraction.sh
