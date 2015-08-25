Compare two BAM files holding alignments of the same reads to different genomes, assigning the alignments to each genome (or none) according to their alignment score. Only single-end files are supported, since this is meant to be used to follow-up fastq_screen results. There is very little error checking performed and only the MAPQ is altered.

The alignments *must* be made with bowtie2. The BAM files *must* have alignments in the same order.

Note that the output files are sorted, so you can directly index and load them in IGV.

Compilation
===========

    Make

Usage
=====

    ./cmpAlignments file1.bam file2.bam

This will write alignments to `file1.unique.bam` and `file2.unique.bam`. The MAPQ scores will be recalculated.
