CC = gcc
OPTS = -Wall -g

.PHONY: all clean htslib clean-all run

.SUFFIXES: .c .o

all: BAM2BedGraph

OBJS = queue.o worker.o

.c.o:
	$(CC) -c $(OPTS) -I../htslib $< -o $@

htslib:
	$(MAKE) -C ../htslib

BAM2BedGraph: htslib $(OBJS)
	$(CC) $(OPTS) -I../htslib -o $@ main.c $(OBJS) ../htslib/libhts.a -lm -lz -lpthread

run: BAM2BedGraph
	time ./BAM2BedGraph "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhRnaSeq/wgEncodeSydhRnaSeqK562Ifna6hPolyaAln.bam" | gzip > foo.cnts.gz

clean:
	rm -f *.o BAM2BedGraph

clean-all: clean
	$(MAKE) --directory=../htslib clean
