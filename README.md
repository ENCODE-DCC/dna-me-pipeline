
DNA methylation pipeline
========================

Before running the pipeline, please edit 'software.path'
to add *absolute* software paths of Bismark and Bowtie.

If you want to leave all outputs from Bismark, please
comment out the last part of the program.


Usage: ./dna-me.sh [options] -i reads.fq -g reference.fa -o outputDir

Required input
--------------
    -i   Single-read fastq input file
    -g   Reference genome sequences
    -o   Output directory

Option
------
    -d   Genome index

    If the above option is specified, the genome indexing step with
    bismark_genome_preparation will be skipped.

Cosmetic Options
----------------
    -h   Show this message
    -q   Don't print what dna-me.sh is doing

Output
------
    Output files will be stored in output directory specified.
    Estimated methylation ratios are in BED files, whose prefixes
    are 'CG_', 'CHG_', and 'CHH_'.

    Each column is composed of:
        [1] chromosome
        [2] start (0-based)
        [3] end   (1-based)
        [4] trinucleotide at mC context
        [5] mC ratio
        [6] strand
