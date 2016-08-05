//gcc -o ~/bin/compareBAMs_withBAMoutput compareBAMs_withBAMoutput.c -I/home/ryand/include -L/home/ryand/lib -lhts -lz -lpthread -Wall
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include "htslib/sam.h"

void usage(char *prog) {
    printf("Usage: %s file1.bam file2.bam ...\n", prog);
    printf("Compare alignments to two organisms, new BAM files containing alignments best assigned to each organism.\n");
}

int main(int argc, char *argv[]) {
    htsFile *fp1, *fp2, *of1, *of2;
    bam_hdr_t *hdr1, *hdr2;
    bam1_t *r1_1 = bam_init1();
    bam1_t *r2_1 = bam_init1();
    bam1_t *r1_2 = bam_init1();
    bam1_t *r2_2 = bam_init1();
    char foo[1024];
    uint8_t *AS1_1, *AS1_2, *AS2_1 = NULL, *AS2_2 = NULL;
    int32_t AS1, AS2;

    if(argc != 3) {
        usage(argv[0]);
        return 1;
    }

    fp1 = sam_open(argv[1], "rb");
    fp2 = sam_open(argv[2], "rb");
    hdr1 = sam_hdr_read(fp1);
    hdr2 = sam_hdr_read(fp2);
    snprintf(foo, 1024, "%s.unique.bam", argv[1]);
    of1 = sam_open(foo, "wbh");
    assert(sam_hdr_write(of1, hdr1) == 0);
    snprintf(foo, 1024, "%s.unique.bam", argv[2]);
    of2 = sam_open(foo, "wbh");
    assert(sam_hdr_write(of2, hdr2) == 0);

    while(sam_read1(fp1, hdr1, r1_1) > 1) {
        if(sam_read1(fp2, hdr2, r1_2) <= 1) break;
        AS1_1 = bam_aux_get(r1_1, "AS");
        AS1_2 = bam_aux_get(r1_2, "AS");

        if(r1_1->core.flag & 1) {
            if(sam_read1(fp1, hdr1, r2_1) < 1) break;
            if(sam_read1(fp2, hdr2, r2_2) < 1) break;
            AS2_1 = bam_aux_get(r2_1, "AS");
            AS2_2 = bam_aux_get(r2_2, "AS");
        }

        //If both pairs can map, then that's preferred
        if((AS1_1 && AS2_1) && (AS1_2 && AS2_2)) {
            AS1 = *((uint8_t*) AS1_1) + *((uint8_t*) AS2_1);
            AS2 = *((uint8_t*) AS1_2) + *((uint8_t*) AS2_2);
            if(AS1 > AS2) {
                assert(sam_write1(of1, hdr1, r1_1));
                assert(sam_write1(of1, hdr1, r2_1));
            } else if(AS2 > AS1) {
                assert(sam_write1(of2, hdr2, r1_2));
                assert(sam_write1(of2, hdr2, r2_2));
            }
        } else if(AS1_1 && AS2_1) {
            assert(sam_write1(of1, hdr1, r1_1));
            assert(sam_write1(of1, hdr1, r2_1));
        } else if(AS1_2 && AS2_2) {
            assert(sam_write1(of2, hdr2, r1_2));
            assert(sam_write1(of2, hdr2, r2_2));
        } else if(AS1_1 && AS1_2) {
            AS1 = *((uint8_t*) AS1_1);
            AS2 = *((uint8_t*) AS1_2);
            if(AS1 > AS2) {
                assert(sam_write1(of1, hdr1, r1_1));
            } else if(AS2 > AS1) {
                assert(sam_write1(of2, hdr2, r1_2));
            }
        } else if(AS1_1) {
            assert(sam_write1(of1, hdr1, r1_1));
        } else if(AS1_2) {
            assert(sam_write1(of2, hdr2, r1_2));
        }
    }

    bam_hdr_destroy(hdr1);
    bam_hdr_destroy(hdr2);
    sam_close(fp1);
    sam_close(fp2);
    sam_close(of1);
    sam_close(of2);

    bam_destroy1(r1_1);
    bam_destroy1(r2_1);
    bam_destroy1(r1_2);
    bam_destroy1(r2_2);
    return 0;
}
