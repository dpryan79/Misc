include=-I/home/ryan/Documents/SourceCode/libGTF -I/home/ryan/Documents/SourceCode/htslib
libs=/home/ryan/Documents/SourceCode/libGTF/libGTF.a /home/ryan/Documents/SourceCode/htslib/libhts.a

.PHONY: all

all: countRepeats

countRepeats:
	gcc -g $(include) -o $@ countRepeats.c $(libs) -lz -lpthread -lpcre

clean:
	rm -f countRepeats
