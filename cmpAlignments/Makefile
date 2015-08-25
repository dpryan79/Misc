CC = gcc
OPTS = -Wall -g

.PHONY: all clean htslib

.SUFFIXES:.c .o

all: cmpAlignments

htslib: 
	$(MAKE) -C ../htslib

cmpAlignments: htslib
	$(CC) $(OPTS) -I../htslib -o cmpAlignments cmpAlignments.c ../htslib/libhts.a -lz -lpthread -lm

clean:
	rm -f cmpAlignments
