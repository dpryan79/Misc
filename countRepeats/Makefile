include=-I../libGTF -I../htslib
libs=../libGTF/libGTF.a ../htslib/libhts.a

.PHONY: all htslib libGTF

all: countRepeats

htslib:
	$(MAKE) -C ../htslib

libGTF:
	$(MAKE) -C ../libGTF

countRepeats: htslib libGTF
	gcc -g $(include) -o $@ countRepeats.c $(libs) -lz -lpthread -lpcre

clean:
	rm -f countRepeats
