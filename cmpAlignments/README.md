Compare two BAM files holding alignments of the same reads to different genomes, assigning the alignments to each genome (or none) according to their alignment score. An inter-file MAPQ threshold can be specified and if the inter-file MAPQ is not at least it then the alignments are ignored. Both paired-end and single-end experiments are supported. There is very little error checking performed.

The alignments *must* be made with bowtie2. The BAM files *must* have alignments in the same order.

Compilation
===========

    Make

Usage
=====

    ./cmpAlignments file1.bam file2.bam

This will write alignments to `file1.unique.bam` and `file2.unique.bam`. To specify an inter-file MAPQ threshold (the default is 0), use the `-q` option:

    ./cmpAlignments -q 5 file1.bam file2.bam
