//gcc -o ~/bin/compareBAMs compareBAMs.c -I/home/ryand/include -L/home/ryand/lib -lhts -lz -lpthread -Wall
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include "htslib/sam.h"

void usage(char *prog) {
    printf("Usage: %s file1.bam file2.bam ...\n", prog);
    printf("Compare alignments to two organisms, outputng the number of alignments to each.\n\
The input is two sets of files.\n");
}

int main(int argc, char *argv[]) {
    htsFile *fp1, *fp2;
    bam_hdr_t *hdr1, *hdr2;
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();
    int i, j, k;
    uint32_t n1, n2;
    char *p;
    uint8_t *aux1, *aux2;

    if(argc < 3) {
        usage(argv[0]);
        return 1;
    }

    for(j=1, k=argc/2+1; k<argc; j++, k++) {
        p = strchr(argv[j], ',');
        if(p) *p = '\0';
        p = strchr(argv[k], ',');
        if(p) *p = '\0';
        fprintf(stderr, "%s vs. %s\n", argv[j], argv[k]);
        n1 = n2 = 0;

        fp1 = sam_open(argv[j], "rb");
        fp2 = sam_open(argv[k], "rb");
        hdr1 = sam_hdr_read(fp1);
        hdr2 = sam_hdr_read(fp2);

        while(sam_read1(fp1, hdr1, read1) > 1) {
            if(sam_read1(fp2, hdr2, read2) <= 1) break;
            aux1 = bam_aux_get(read1, "AS");
            aux2 = bam_aux_get(read2, "AS");
            if(aux1 && !aux2) {
                if(read1->core.qual < 2) continue;
                n1++;
             } else if(aux2 && !aux1) {
                if(read2->core.qual < 2) continue;
                n2++;
             }
        }

        bam_hdr_destroy(hdr1);
        bam_hdr_destroy(hdr2);
        sam_close(fp1);
        sam_close(fp2);

        p = strchr(argv[j], '.');
        if(p) *p = '\0';
        printf("%s\t%"PRIu32"\t%"PRIu32"\n", argv[j], n1, n2);
    }

    bam_destroy1(read1);
    bam_destroy1(read2);
    return 0;
}
