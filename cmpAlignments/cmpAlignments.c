#include "htslib/sam.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

int32_t getAS(bam1_t *b) {
    int32_t def = -1073741823; //half of INT_MIN
    uint8_t *p = bam_aux_get(b, "AS");
    if(!p) return def;
    return bam_aux2i(p);
}

int32_t getXS(bam1_t *b) {
    int32_t def = -1073741823;
    uint8_t *p = bam_aux_get(b, "XS");
    if(!p) return def;
    return bam_aux2i(p);
}

int32_t calcScMin(bam1_t *b) {
    return -0.6-0.6*b->core.l_qseq;
}

//This is just copied from Bison
int32_t calcMAPQ(int32_t AS, int32_t XS, int32_t scMin) {
    int diff, bestOver, bestdiff;
    diff = abs(scMin); //Range of possible alignment scores
    bestOver = AS-scMin; //Shift alignment score range, so worst score is 0
    
    bestdiff = (int) abs(abs((float) AS)-abs((float) XS)); //Absolute distance between alignment scores
    if(XS < scMin) {
        if(bestOver >= diff * (double) 0.8f) return 42;
        else if(bestOver >= diff * (double) 0.7f) return 40;
        else if(bestOver >= diff * (double) 0.6f) return 24;
        else if(bestOver >= diff * (double) 0.5f) return 23;
        else if(bestOver >= diff * (double) 0.4f) return 8;
        else if(bestOver >= diff * (double) 0.3f) return 3;
        else return 0;
    } else {
        if(bestdiff >= diff * (double) 0.9f) {
            if(bestOver == diff) {
                return 39;
            } else {
                return 33;
            }
        } else if(bestdiff >= diff * (double) 0.8f) {
            if(bestOver == diff) {
                return 38;
            } else {
                return 27;
            }
        } else if(bestdiff >= diff * (double) 0.7f) {
            if(bestOver == diff) {
                return 37;
            } else {
                return 26;
            }
        } else if(bestdiff >= diff * (double) 0.6f) {
            if(bestOver == diff) {
                return 36;
            } else {
                return 22;
            }
        } else if(bestdiff >= diff * (double) 0.5f) {
            if(bestOver == diff) {
                return 35;
            } else if(bestOver >= diff * (double) 0.84f) {
                return 25;
            } else if(bestOver >= diff * (double) 0.68f) {
                return 16;
            } else {
                return 5;
            }
        } else if(bestdiff >= diff * (double) 0.4f) {
            if(bestOver == diff) {
                return 34;
            } else if(bestOver >= diff * (double) 0.84f) {
                return 21;
            } else if(bestOver >= diff * (double) 0.68f) {
                return 14;
            } else {
                return 4;
            }
        } else if(bestdiff >= diff * (double) 0.3f) {
            if(bestOver == diff) {
                return 32;
            } else if(bestOver >= diff * (double) 0.88f) {
                return 18;
            } else if(bestOver >= diff * (double) 0.67f) {
                return 15;
            } else {
                return 3;
            }
        } else if(bestdiff >= diff * (double) 0.2f) {
            if(bestOver == diff) {
                return 31;
            } else if(bestOver >= diff * (double) 0.88f) {
                return 17;
            } else if(bestOver >= diff * (double) 0.67f) {
                return 11;
            } else {
                return 0;
            }
        } else if(bestdiff >= diff * (double) 0.1f) {
            if(bestOver == diff) {
                return 30;
            } else if(bestOver >= diff * (double) 0.88f) {
                return 12;
            } else if(bestOver >= diff * (double) 0.67f) {
                return 7;
            } else {
                return 0;
            }
        } else if(bestdiff > 0) {
            if(bestOver >= diff * (double)0.67f) {
                return 6;
            } else {
                return 2;
            }
        } else {
            if(bestOver >= diff * (double)0.67f) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}

int writeAlignment(htsFile *o, bam_hdr_t *hdr, bam1_t *b, int32_t MAPQ) {
    assert(MAPQ >= 0);
    assert(MAPQ <= 255);
    b->core.qual = MAPQ;
    return sam_write1(o, hdr, b);
}

int main(int argc, char *argv[]) {
    htsFile *i1, *i2, *o1, *o2, *o;
    bam1_t *r1, *r2;
    bam_hdr_t *hdr1, *hdr2;
    int32_t AS1, AS2, XS1, XS2;
    int32_t bestAS, bestXS, MAPQ;
    int32_t scMin = -1073741823;
    char str[1024], *p;

    if(argc != 3) {
        fprintf(stderr, "Usage: %s file1.bam file2.bam\n", argv[0]);
        fprintf(stderr, "Unique alignments with recalculated MAPQ scores will be written to file1.unique.bam and file2.unique.bam\n");
        return 0;
    }

    i1 = sam_open(argv[1], "rb");
    assert(i1);
    i2 = sam_open(argv[2], "rb");
    assert(i2);

    hdr1 = sam_hdr_read(i1);
    assert(hdr1);
    hdr2 = sam_hdr_read(i2);
    assert(hdr2);

    r1 = bam_init1();
    r2 = bam_init1();
    assert(r1);
    assert(r2);

    //create the output file
    p = strrchr(argv[1], '.');
    if(p) p[0] = '\0';
    p = strrchr(argv[2], '.');
    if(p) p[0] = '\0';
    snprintf(str, 1024, "%s.unique.bam", argv[1]);
    o1 = sam_open(str, "wb");
    snprintf(str, 1024, "%s.unique.bam", argv[2]);
    o2 = sam_open(str, "wb");

    assert(sam_hdr_write(o1, hdr1));
    assert(sam_hdr_write(o2, hdr2));

    while(sam_read1(i1, hdr1, r1) > 0) {
        assert(sam_read1(i1, hdr2, r2) > 0);

        //Make sure the names match
        assert(strcmp(bam_get_qname(r1), bam_get_qname(r2)) == 0);

        //Get the scores
        if(scMin == -1073741823) scMin = calcScMin(r1);
        AS1 = getAS(r1);
        AS2 = getAS(r2);
        XS1 = getXS(r1);
        XS2 = getXS(r2);

        if(AS1 > AS2) {//Assign to o1
            bestAS = AS1;
            if(XS1 > AS2) bestXS = XS1;
            else bestXS = AS2;
            o = o1;
        } else if(AS2 > AS1) { //Assign to o2
            bestAS = AS2;
            if(AS1 > XS2) bestXS = AS1;
            else bestXS = XS2;
            o = o2;
        } else {
            continue;
        }

        MAPQ = calcMAPQ(bestAS, bestXS, scMin);
        if(o == o1) {
            assert((r1->core.flag&4) == 0);
            assert(writeAlignment(o1, hdr1, r1, MAPQ));
        } else {
            assert((r2->core.flag&4) == 0);
            assert(writeAlignment(o2, hdr2, r2, MAPQ));
        }
    }

    bam_destroy1(r1);
    bam_destroy1(r2);
    bam_hdr_destroy(hdr1);
    bam_hdr_destroy(hdr2);
    sam_close(o1);
    sam_close(o2);
    sam_close(i1);
    sam_close(i2);
    return 0;
}
