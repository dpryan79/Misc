This is intended to test an issue in deeptools when remote files are being used. Currently, `bamCoverage` will issue a whole slew of errors and never complete if a remote BAM file on an http server is being used. The errors are printed from pysam (and therefore htslib in [knetfile.c](https://github.com/dpryan79/htslib/blob/develop/knetfile.c#L561)). This only occurs when multiple subprocesses are used via the pool of workers in the "multiprocessing" package.

The question, then is whether this is:
 1. An innate issue with processing BAM files over http
 2. A strange interaction between using BAM files over https and the "multiprocessing" module.

Making a pthreads version in C should clarify this completely.

Results: The warning message is still printed but everything completes as expected, though slowly. There's no reason to use multiple threads for this task, since network IO is too slow for even a single core to be fully used in most cases. Consequently, the deeptools mapReduce paradigm can be dispensed with. The warning message will still get printed, but only once. I can easily put a pull request in to have that `fprintf` removed, since it's correct that you can't seek to the end (I tried).
