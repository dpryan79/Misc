`fastqIntegrity` is aimed at ensuring that an untrimmed fastq file is more or less properly formatted. We recently had an issue where an L3 cache for a CPU on the server where our demultiplexing is done. The following is done for each fastq file:

 * Ensures that each read name begins with '@'
 * Ensures that all base calls are legitimate (i.e., one of A, C, G, T, or N).
 * Ensures that all reads have the same length
 * Ensures that the phred scores are all valid for Illumina 1.8+ format

Note that you can pass in multiple fastq files and then the program will ensure that *all of them have identical read lengths*.

Note also that this currently supports reads shorter than 1023 bases only (simply for the sake of simplicity).

Compilation
===========

    make

Note that this will run tests on each of the 7 test files. Only one of these is correctly formatted and the each test incorrectly formatted fastq files.

Output
======

If the program sees an error, it will return one of the following codes:

    -1: Gzip error
    -2: Name doesn't start with @
    -3: Sequence has an illegal character
    -4: Line 3 doesn't start with '+'
    -5: QUAL has unexpected length
    -6: Unexpected QUAL score
    -7: Incorrect number of lines

If a file is OK, it will report the number of reads encountered, the read length and the file name.
