CC = gcc
OPTS = -Wall -g

.PHONY: all clean htslib clean-all

.SUFFIXES: .c .o

all: fastqIntegrity

.c.o:
	$(CC) -c $(OPTS) -I../htslib $< -o $@

htslib:
	$(MAKE) -C ../htslib

fastqIntegrity: htslib
	$(CC) $(OPTS) -I../htslib -o $@ fastqIntegrity.c ../htslib/libhts.a -lm -lz -lpthread

clean:
	rm -f *.o fastqIntegrity

clean-all: clean
	$(MAKE) --directory=../htslib clean
