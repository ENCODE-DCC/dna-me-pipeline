
DNA methylation pipeline
========================

Before running the pipeline, please edit 'software.path'
to add *absolute* software paths of Bismark and Bowtie.

If you want to leave all outputs from Bismark, please
comment out the last part of the program.


Usage: ./dna-me.sh [options] -g reference.fa -o outputDir {-i <single.fq> | -1 <pair1.fq> -2 <pair2.fq>}

Required input
--------------
    -i   Single-read fastq input file
    -1   Paired-end fastq file containing the #1 mates
    -2   Paired-end fastq file containing the #2 mates
    -g   Reference genome sequences
    -o   Output directory

Option
------
    -d   Genome index

    If the above option is specified, the genome indexing step with
    bismark_genome_preparation will be skipped.

    -I   Minimum insert size for paired-end alignments (default=0)
    -X   Maximum insert size for paired-end alignments (default=500)

    From the Bismark man page (e.g. '-I' option):
      If -I 60 is specified and a paired-end alignment consists of two
      20-bp alignments in the appropriate orientation with a 20-bp gap
      between them, that alignment is considered valid (as long as -X
      is also satisfied). A 19-bp gap would not be valid in that case.

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
        [7] mC reads
        [8] non-mC reads (i.e. all reads - mC reads)
