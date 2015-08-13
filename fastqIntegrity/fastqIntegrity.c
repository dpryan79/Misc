#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>
#include <inttypes.h>

#define BUFLEN 1024

/* returns:
     0: EOF
     1: Continue
    -1: Gzip error
    -2: Name doesn't start with @
    -3: Sequence has an illegal character
    -4: Line 3 doesn't start with '+'
    -5: QUAL has unexpected length
    -6: Unexpected QUAL score
*/
int getRead(char *buf, gzFile fp, int *expectedLength) {
    int i, err;
    //Read name
    buf = gzgets(fp, buf, BUFLEN);
    gzerror(fp, &err);
    if(err != Z_OK) {
        return -1;
    } else if(gzeof(fp) == 1) {
        return 0;
    }
    if(buf[0] != '@') return -2;
    //Sequence
    buf = gzgets(fp, buf, BUFLEN);
    gzerror(fp, &err);
    if(err != Z_OK) {
        return -1;
    } else if(gzeof(fp) == 1) {
        return -7;
    }
    if(*expectedLength) assert(strlen(buf)-1 == *expectedLength);
    else *expectedLength = strlen(buf)-1;
    for(i=0; i<*expectedLength; i++) {
        if((buf[i] != 'A') && (buf[i] != 'C') && (buf[i] != 'G') && (buf[i] != 'T') && (buf[i] != 'N')) return -3;
    }
    //+
    buf = gzgets(fp, buf, BUFLEN);
    gzerror(fp, &err);
    if(err != Z_OK) {
        return -1;
    } else if(gzeof(fp) == 1) {
        return -7;
    }
    if(buf[0] != '+') return -4;
    //Qual
    buf = gzgets(fp, buf, BUFLEN);
    gzerror(fp, &err);
    if(err != Z_OK) {
        return -1;
    } else if(gzeof(fp) == 1) {
        return -7;
    }
    if(strlen(buf)-1 != *expectedLength) return -5;
    for(i=0; i<*expectedLength; i++) {
        if((buf[i] < 33) || (buf[i] > 73)) return -6;
    }

    return 1;
}

int main(int argc, char *argv[]) {
    gzFile fp;
    char buf[BUFLEN];
    int expectedLength = 0, i, rv;
    uint32_t total = 0;

    if(argc<2) {
        fprintf(stderr, "Usage: %s file.fastq.gz\n", argv[0]);
        return 1;
    }

    for(i=1; i<argc; i++) {
        fp = gzopen(argv[i], "r");
        if(!fp) {
            fprintf(stderr, "Couldn't open %s for reading!\n", argv[1]);
            return 1;
        }

        total = 0;
        while((rv = getRead(buf, fp, &expectedLength)) == 1) total++;
        if(rv != 0) printf("An error occured (%i) while processing %s!\n", rv, argv[i]);
        else printf("%"PRIu32" reads of length %i %s\n", total, expectedLength, argv[i]);

        gzclose(fp);
    }

    return 0;
}
