prefix= /usr/local/bin
CC = gcc
OPTS = -Wall -g

.PHONY: all clean htslib install clean-all

all: countInternalBarcodes

htslib: 
	$(MAKE) -C ../htslib

countInternalBarcodes: htslib
	$(CC) $(OPTS) -I../htslib -o countInternalBarcodes countInternalBarcodes.c ../htslib/libhts.a -lm -lz -lpthread

clean:
	rm -f countInternalBarcodes

clean-all: clean
	make --directory=../htslib clean

install: countInternalBarcodes
	install countInternalBarcodes $(prefix)
