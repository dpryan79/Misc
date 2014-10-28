//gcc -o ~/bin/CompareBAMs CompareBAMs.c -I/home/ryand/include -L/home/ryand/lib -lhts -lz -lpthread -Wall
#include <stdio.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"

void usage(char *prog) {
    printf("Usage: %s file1.bam file2.bam\n", prog);
    printf("\n\
Compare two BAM files, alignment by alignment. Only alignments with identical\n\
names will be compared. Differences can occur due to (1) tid, (2) position, and\n\
(3) CIGAR operations. This is useful when comparing pre- and post-realignment.\n");
}

int main(int argc, char *argv[]) {
    htsFile *fp1, *fp2;
    bam_hdr_t *hdr1, *hdr2;
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();
    kstring_t *ks = calloc(1, sizeof(kstring_t));
    int i;
    uint32_t *cigar1, *cigar2;

    if(argc != 3) {
        usage(argv[0]);
        return 1;
    }
    fp1 = sam_open(argv[1], "rb");
    fp2 = sam_open(argv[2], "rb");
    hdr1 = sam_hdr_read(fp1);
    hdr2 = sam_hdr_read(fp2);

    while(sam_read1(fp1, hdr1, read1) > 1) {
        if(sam_read1(fp2, hdr2, read2) <= 1) break;
        if(strcmp(bam_get_qname(read1), bam_get_qname(read2)) == 0) {
            if(read1->core.tid != read2->core.tid) {
                printf("Differ on tid!\n");
                sam_format1(hdr1, read1, ks);
                printf("%s\n", ks->s);
                sam_format1(hdr2, read2, ks);
                printf("%s\n", ks->s);
            } else {
                //Start
                if(read1->core.pos != read2->core.pos) {
                    printf("Differ on pos!\n");
                    sam_format1(hdr1, read1, ks);
                    printf("%s\n", ks->s);
                    sam_format1(hdr2, read2, ks);
                    printf("%s\n", ks->s);
                } else {
                    if(read1->core.n_cigar != read2->core.n_cigar) {
                        printf("Differ on n_cigar!\n");
                        sam_format1(hdr1, read1, ks);
                        printf("%s\n", ks->s);
                        sam_format1(hdr2, read2, ks);
                        printf("%s\n", ks->s);
                    } else {
                        cigar1 = bam_get_cigar(read1);
                        cigar2 = bam_get_cigar(read2);
                        for(i=0; i<read1->core.n_cigar; i++) {
                            if(cigar1[i] != cigar2[i]) {
                                printf("Differ on cigar operation %i!\n", i);
                                sam_format1(hdr1, read1, ks);
                                printf("%s\n", ks->s);
                                sam_format1(hdr2, read2, ks);
                                printf("%s\n", ks->s);
                                break;
                            }
                        } //cigar ops
                    } //n_cigar
                } //pos
            } //tid
        } //qname
    }

    bam_hdr_destroy(hdr1);
    bam_hdr_destroy(hdr2);
    sam_close(fp1);
    sam_close(fp2);
    bam_destroy1(read1);
    bam_destroy1(read2);
    free(ks);
    return 0;
}
