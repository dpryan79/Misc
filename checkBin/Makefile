CC = gcc
OPTS = -Wall -g

.PHONY: all clean htslib clean-all

.SUFFIXES: .c .o

all: checkBin

.c.o:
	$(CC) -c $(OPTS) -I../htslib $< -o $@

htslib:
	$(MAKE) -C ../htslib

checkBin: htslib
	$(CC) $(OPTS) -I../htslib -o $@ checkBin.c ../htslib/libhts.a -lm -lz -lpthread

clean:
	rm -f *.o checkBin

clean-all: clean
	$(MAKE) --directory=../htslib clean
