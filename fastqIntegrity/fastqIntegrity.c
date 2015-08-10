#include "../htslib/htslib/kseq.h"
#include "../htslib/htslib/kstring.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>
#include <inttypes.h>

KSTREAM_INIT(gzFile, gzread, 16384)

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
int getRead(kstring_t *ks, kstream_t *seq, int *expectedLength) {
    int ret, i, err;
    //Read name
    ks_getuntil2(seq, '\n', ks, &ret, 0);
    if(ret == 0) {
        if(gzeof(seq->f) == 1) {
            return 0;
        } else {
            fprintf(stderr, "Got gzerror %s\n", gzerror(seq->f, &err));
            return -1;
        }
    }
    if(ks->s[0] != '@') return -2;
    //Sequence
    ks_getuntil2(seq, '\n', ks, &ret, 0);
    if(ret == 0) {
        if(gzeof(seq->f) == 1) {
            return 0;
        } else {
            fprintf(stderr, "Got gzerror %s\n", gzerror(seq->f, &err));
            return -1;
        }
    }
    if(*expectedLength) assert(strlen(ks->s) == *expectedLength);
    else *expectedLength = strlen(ks->s);
    for(i=0; i<*expectedLength; i++) {
        if((ks->s[i] != 'A') && (ks->s[i] != 'C') && (ks->s[i] != 'G') && (ks->s[i] != 'T') && (ks->s[i] != 'N')) return -3;
    }
    //+
    ks_getuntil2(seq, '\n', ks, &ret, 0);
    if(ret == 0) {
        if(gzeof(seq->f) == 1) {
            return 0;
        } else {
            fprintf(stderr, "Got gzerror %s\n", gzerror(seq->f, &err));
            return -1;
        }
    }
    if(ks->s[0] != '+') return -4;
    //Qual
    ks_getuntil2(seq, '\n', ks, &ret, 0);
    if(ret == 0) {
        if(gzeof(seq->f) == 1) {
            return 0;
        } else {
            fprintf(stderr, "Got gzerror %s\n", gzerror(seq->f, &err));
            return -1;
        }
    }
    if(strlen(ks->s) != *expectedLength) return -5;
    for(i=0; i<*expectedLength; i++) {
        if((ks->s[i] < 33) || (ks->s[i] > 73)) return -6;
    }

    return 1;
}

int main(int argc, char *argv[]) {
    gzFile fp;
    kstream_t *kstream;
    kstring_t *ks;
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

        kstream = ks_init(fp);
        assert(kstream);
        ks = calloc(1, sizeof(kstring_t));
        assert(ks);
        total = 0;
        while((rv = getRead(ks, kstream, &expectedLength)) == 1) total++;
        if(rv != 0) printf("An error occured (%i) while processing %s!\n", rv, argv[i]);
        else printf("%"PRIu32" reads of length %i %s\n", total, expectedLength, argv[i]);

        ks_destroy(kstream);
        free(ks->s);
        free(ks);
        gzclose(fp);
    }

    return 0;
}
