#include "htslib/sam.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>

int32_t getAS(bam1_t *b) {
    int32_t def = -1073741823; //half of INT_MIN
    uint8_t *p = bam_aux_get(b, "AS");
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

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] file1.bam file2.bam\n", prog);
    printf("Alignments in each file are compared, line by line, and those pairs with an\n"
           "INTER-file MAPQ above threshold are written to file1.unique.bam or\n"
           "file2.unique.bam, as appropriate.\n"
           "\nOptions:\n"
           " -q INT  The inter-file MAPQ threshold to write a pair to the output. Default 0.\n"
           " -h      Print this help message.\n");
}

int main(int argc, char *argv[]) {
    htsFile *i1, *i2, *o1, *o2, *o;
    bam1_t *r1_1, *r1_2, *r2_1, *r2_2;
    bam_hdr_t *hdr1, *hdr2;
    int32_t AS1_1, AS1_2, AS2_1, AS2_2, AS, XS;
    int32_t scMin1 = -1073741823;
    int32_t scMin2 = -1073741823;
    char str[1024], *p;
    int opt, minMAPQ = 0;

    while((opt = getopt(argc, argv, "hq:")) != -1) {
        switch(opt) {
        case 'h':
            usage(argv[0]);
            return 0;
        case 'q':
            minMAPQ = atoi(optarg);
            break;
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if(argc-optind != 2) {
        usage(argv[0]);
        return 0;
    }

    i1 = sam_open(argv[optind], "r");
    assert(i1);
    i2 = sam_open(argv[optind+1], "r");
    assert(i2);

    hdr1 = sam_hdr_read(i1);
    assert(hdr1);
    hdr2 = sam_hdr_read(i2);
    assert(hdr2);

    r1_1 = bam_init1();
    r1_2 = bam_init1();
    r2_1 = bam_init1();
    r2_2 = bam_init1();
    assert(r1_1);
    assert(r1_2);
    assert(r2_1);
    assert(r2_2);

    //create the output file
    p = strrchr(argv[optind], '.');
    if(p) p[0] = '\0';
    p = strrchr(argv[optind+1], '.');
    if(p) p[0] = '\0';
    snprintf(str, 1024, "%s.unique.bam", argv[optind]);
    o1 = sam_open(str, "wb");
    snprintf(str, 1024, "%s.unique.bam", argv[optind+1]);
    o2 = sam_open(str, "wb");

    assert(sam_hdr_write(o1, hdr1) == 0);
    assert(sam_hdr_write(o2, hdr2) == 0);

    while(sam_read1(i1, hdr1, r1_1) > 0) {
        if(r1_1->core.flag & 1) assert(sam_read1(i1, hdr1, r1_2));

        assert(sam_read1(i2, hdr2, r2_1) > 0);
        if(r2_1->core.flag & 1) assert(sam_read1(i2, hdr2, r2_2));

        //Make sure the names match
        assert(strcmp(bam_get_qname(r1_1), bam_get_qname(r2_1)) == 0);
        if(r1_1->core.flag & 1) {
            assert(strcmp(bam_get_qname(r1_1), bam_get_qname(r1_2)) == 0);
            assert(strcmp(bam_get_qname(r1_1), bam_get_qname(r2_2)) == 0);
        }

        //Get the scores
        scMin1 = calcScMin(r1_1);
        scMin2 = 0;
        if(r1_1->core.flag & 1) scMin2 = calcScMin(r1_2);
        AS1_1 = getAS(r1_1);
        AS2_1 = getAS(r2_1);
        if(r1_1->core.flag & 1) {
            AS1_2 = getAS(r1_2);
            AS2_2 = getAS(r2_2);
        } else {
            AS1_2 = AS2_2 = 0;
        }

        if(AS1_1+AS1_2 > AS2_1+AS2_2) {//Assign to o1
            AS = AS1_1 + AS1_2;
            XS = AS2_1 + AS2_2;
            o = o1;
        } else if(AS2_1+AS2_2 > AS1_1+AS1_2) { //Assign to o2
            AS = AS2_1 + AS2_2;
            XS = AS1_1 + AS1_2;
            o = o2;
        } else {
            continue;
        }

        if(calcMAPQ(AS, XS, scMin1+scMin2) >= minMAPQ) {
            if(o == o1) {
                assert(sam_write1(o, hdr1, r1_1));
                if(r1_1->core.flag & 1) assert(sam_write1(o, hdr1, r1_2));
            } else {
                assert(sam_write1(o, hdr2, r2_1));
                if(r2_1->core.flag & 1) assert(sam_write1(o, hdr2, r2_2));
            }
        }
    }

    bam_destroy1(r1_1);
    bam_destroy1(r1_2);
    bam_destroy1(r2_1);
    bam_destroy1(r2_2);
    bam_hdr_destroy(hdr1);
    bam_hdr_destroy(hdr2);
    sam_close(o1);
    sam_close(o2);
    sam_close(i1);
    sam_close(i2);
    return 0;
}
