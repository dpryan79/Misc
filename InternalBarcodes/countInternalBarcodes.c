#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include <getopt.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/kseq.h"

//For convenience
int nBarcodes=36;
char *barcodes = "ATCACGCGATGTTTAGGCTGACCAACAGTGGCCAATCAGATCACTTGAGATCAGTAGCTTGGCTACCTTGTAAGTCAAAGTTCCATGTCACCGTCCGTCCGCGTGAAAGTGGCCGTTTCGCGTACGGAGTGGACTGATATTCCTCGTGATACATCGGCCTAATGGTCACACTGTATTGGCGATCTGTCAAGTCTGATCAAGCTAGTAGCCTACAAG";

KSTREAM_INIT2(, gzFile, gzread, 16384)

typedef struct {
    int64_t score;
    int32_t end;
    int adapterType; //0: Truseq; 1: NEBNext
} s_align;

/*******************************************************************************
*
* This is modified from MethIndelRealigner
*
* The memory demands could be lowered to 4*refLen
*
*******************************************************************************/
s_align * GlobalAlignment(char *ref, int32_t refLen, char *seq, int32_t seqLen) {
    int64_t i, j, *last, *cur, *tmp;
    int8_t **mat;
    int64_t nmatch = -1, mismatch = -3, gapOpen = -5, gapExtend = -3, left, top, diag, best;
    s_align *sal = malloc(sizeof(s_align));
    assert(sal);

    sal->score = -3*refLen;
    sal->end = -1;

    //Initialize the score vectors
    mat = malloc((seqLen+1)*sizeof(int8_t*));
    assert(mat);
    for(i=0; i<seqLen+1; i++) {
        mat[i] = calloc(refLen+1, sizeof(int8_t));
        assert(mat[i]);
    }
    for(j=1; j<refLen+1; j++) mat[0][j] = 4; //left
    last = calloc(refLen+1, sizeof(int64_t));
    cur = malloc(sizeof(int64_t) * (refLen+1));
    assert(last); assert(cur);

    //Fill in the direction matrix
    for(i=1; i<seqLen+1; i++) {
        mat[i][0] = 1;
        cur[0] = gapExtend*i + gapOpen;
        for(j=1; j<refLen+1; j++) {
            left = cur[j-1] + gapExtend + ((mat[i][j-1] & 4)?0:gapOpen);
            top = last[j] + gapExtend + ((mat[i-1][j] & 1)?0:gapOpen);
            diag = last[j-1];
            if(seq[i-1] != ref[j-1]) {
                if(seq[i-1] == 'N') diag += nmatch;
                else diag += mismatch;
            }
            best = left;
            if(best < top) best = top;
            if(best < diag) best = diag;
            cur[j] = best;
            if(best == top) mat[i][j] |= 1;
            if(best == diag) mat[i][j] |= 2;
            if(best == left) mat[i][j] |= 4;
        }
        tmp = cur;
        cur = last;
        last = tmp;

        //If the most positive value is < -12 then break
        left = cur[0];
        for(j=1; j<refLen+1; j++) {
            if(cur[j] > left) left=cur[j];
        }
        if(left < -12) break;
    }

    for(i=1; i<refLen+1; i++) {
        if(cur[i] >= sal->score) { //take the right-most case
            sal->score = cur[i];
            sal->end = i;
        }
    }

    //Aufraumen! -Tristan
    for(i=0; i<seqLen+1; i++) free(mat[i]);
    free(mat);
    free(cur);
    free(last);

    return sal;
}

char RC(char c) {
    switch(c) {
    case 'A' :
        return 'T';
    case 'C' :
        return 'G';
    case 'G' :
        return 'C';
    case 'T' :
        return 'A';
    default :
        return 'N';
    }
}

void reverseComp(kstring_t *read) {
    char c;
    int i;

    for(i=0; i<(read->l>>1) + (read->l&1); i++) {
        c = RC(read->s[i]);
        read->s[i] = RC(read->s[read->l-1-i]);
        read->s[read->l-1-i] = c;
    }
}

//Return the barcode index
//Allow an alignment score of up to -12
int hasBarcode(kstring_t *read) {
    char *adapter = "CAAGCAGAAGACGGCATACGAGAT"; //NEBNext
    char *adapter2 ="AGCACACGTCTGAACTCCAGTCAC"; //TruSeq
    int i, start, end;

    //NEBNext
    //Allow a score of up to -12, i.e., 4 mismatches or a indel plus mismatch
    s_align *s = GlobalAlignment(read->s, read->l, adapter, 24);
    s->adapterType = 1;
    if(s->score < -12) { //TruSeq
        free(s);
        s = GlobalAlignment(read->s, read->l, adapter2, 24);
        s->adapterType = 0;
    }
    if(s->score >= -12) {
        if(s->adapterType == 0) {
            start = 24;
            end = nBarcodes;
        } else if(s->adapterType == 1) {
            start = 0;
            end = 24;
        }
        for(i=start; i<end; i++) {
            if(strncmp(barcodes+6*i, read->s+s->end+1, 6) == 0) {
                free(s);
                return i;
            }
        }
    }
    free(s);

    //What about the reverse complement?
    reverseComp(read);

    //NEBNext
    s = GlobalAlignment(read->s, read->l, adapter, 24);
    s->adapterType = 1;
    if(s->score < -12) { //TruSeq
        free(s);
        s = GlobalAlignment(read->s, read->l, adapter2, 24);
        s->adapterType = 0;
    }
    if(s->score >= -12) {
        if(s->adapterType == 0) {
            start = 24;
            end = nBarcodes;
        } else if(s->adapterType == 1) {
            start = 0;
            end = 24;
        }
        for(i=start; i<end; i++) {
            if(strncmp(barcodes+6*i, read->s+s->end+1, 6) == 0) {
                free(s);
                return i;
            }
        }
    }
    free(s);
    return -1;
}

int bcIDX(kstring_t ks) {
    int i;
    char *p = strstr(ks.s, " 1:N:0:");
    if(p) p += 7;
    for(i=0; i<nBarcodes; i++) if(strncmp(barcodes+6*i, p, 6) == 0) return i;
    return -1;
}

void usage() {
    fprintf(stderr, "Usage: countInternalBarcodes <file.fastq.gz>\n");
}

int main(int argc, char *argv[]) {
    gzFile gz;
    char c;
    int ret, idx, bc = -1;
    kstream_t *ks;
    kstring_t line1, line2;
    uint32_t *cnts;

    line1.s = 0; line1.l = line1.m = 0;
    line2.s = 0; line2.l = line2.m = 0;

    static struct option lopts[] = {
        {"help", 0, NULL, 'h'},
        {0,      0, NULL,   0}
    };
    while((c = getopt_long(argc, argv, "h", lopts, NULL)) >= 0) {
        switch(c) {
        case 'h' :
            usage();
            return 0;
        default :
            fprintf(stderr, "Invalid option '%c'\n", c);
            usage();
            return 1;
        }
    }

    if(argc == 1) {
        usage();
        return 0;
    }
    if(argc-optind != 1) {
        fprintf(stderr, "You must provide a fastq file!\n");
        usage();
        return 1;
    }

    gz = gzopen(argv[optind], "rb");
    assert(gz);
    ks = ks_init(gz);
    assert(ks);

    cnts = calloc(nBarcodes, sizeof(uint32_t));

    while(ks_getuntil(ks, KS_SEP_LINE, &line1, &ret) >= 0) {
        if(bc == -1) bc = bcIDX(line1);
        assert(ks_getuntil(ks, KS_SEP_LINE, &line2, &ret) >= 0);
        idx = hasBarcode(&line2);
        if(idx>=0) cnts[idx] += 1;
        assert(ks_getuntil(ks, KS_SEP_LINE, &line2, &ret) >= 0);
        assert(ks_getuntil(ks, KS_SEP_LINE, &line2, &ret) >= 0);
    }

    if(bc>=0) cnts[bc] = 0;
    printf("%"PRIu32, cnts[0]);
    for(ret=1; ret<nBarcodes; ret++) printf("\t%"PRIu32, cnts[ret]);
    printf("\n");

    ks_destroy(ks);
    free(line1.s);
    free(line2.s);
    free(cnts);
    gzclose(gz);

    return 0;
}
