This is intended to test an issue in deeptools when remote files are being used. Currently, `bamCoverage` will issue a whole slew of errors and never complete if a remote BAM file on an http server is being used. The errors are printed from pysam (and therefore htslib in [knetfile.c](https://github.com/dpryan79/htslib/blob/develop/knetfile.c#L561)). This only occurs when multiple subprocesses are used via the pool of workers in the "multiprocessing" package.

The question, then is whether this is:
 1. An innate issue with processing BAM files over http
 2. A strange interaction between using BAM files over https and the "multiprocessing" module.

Make a pthreads version in C should clarify this completely.
