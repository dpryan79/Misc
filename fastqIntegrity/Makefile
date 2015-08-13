CC = gcc
OPTS = -Wall -g

.PHONY: all clean test

.SUFFIXES: .c .o

all: fastqIntegrity test

.c.o:
	$(CC) -c $(OPTS) $< -o $@

fastqIntegrity:
	$(CC) $(OPTS) -o $@ fastqIntegrity.c ../htslib/libhts.a -lz

test:
	./fastqIntegrity tests/*.gz

clean:
	rm -f *.o fastqIntegrity
