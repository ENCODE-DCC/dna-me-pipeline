[Description]
Basically, this is a bismark_pipeline for paired-end fastqs. The idea would be to split the fastq files in to smaller ones (around 18 million reads) and parallelize the job for QC using trim galore, and then simultaneously align in parallel with the bowtie2 instances on different nodes, which would  provide us the bam alignment quicker. Eventually, we could merge those unsorted bam files back to remove the PCR duplicates, and further run on methylation_extractor or coverage metrics on those merged de-duplicated bams.


[Usage:]
./bismark_pipeline_main.sh 

[List of variables]:
Following variables must be changed in 'bismark_pipeline_main.sh'

OUTPUT_LOC, 
INPUT_LOC, 
LIB_LIST, 
GENOME_PATH, 
BISMARK_PATH, 
SAMTOOLS_PATH, 
TRIMGALORE_PATH, 
BOWTIE_PATH, 
CORE_NUM
