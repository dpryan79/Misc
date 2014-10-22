#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "bam.h"

typedef struct {
    int32_t id;
    int32_t tid;
    int32_t pos1, pos2;
} reads;

inline int cmp_reads(const void *a, const void *b) {
    reads *r1 = (reads *) a;
    reads *r2 = (reads *) b;

    return((r2->id) - (r1->id));
}

inline int compare_reads(reads *r1, reads *r2) {
    int rv = (r2->id) - (r1->id);

    if(rv != 0) return rv; //Different reads
    if(r1->tid == r2->tid && r1->pos1 == r2->pos1 && r1->pos2 == r2->pos2) return 0; //The same mapping
    else return -2; //Different mapping
}

int main(int argc, char *argv[]) {
    unsigned long length1 = 0, length2 = 0;
    unsigned long max_length1 = 1000000, max_length2 = 1000000;
    unsigned long i = 0, j = 0, n1 = 0, n2 = 0;
    bamFile fp1, fp2;
    bam_header_t *header1, *header2;
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();
    reads *reads1 = malloc(max_length1 * sizeof(reads));
    reads *reads2 = malloc(max_length2 * sizeof(reads));
    char *name, *p, *oname;
    FILE *of;
    int rv;

    if(argc != 3) {
        printf("Usage: %s file1.bam file2.bam\n", argv[0]);
        return 0;
    }

    //Create the output name
    oname = malloc(sizeof(char) * (strlen(argv[1]) + strlen(argv[2]) + strlen("_vs_.txt ")));
    sprintf(oname, "%s_vs_%s.txt", argv[1], argv[2]);
    of = fopen(oname, "w");
    printf("Output will be written to %s\n", oname);

    //Open the inputs
    fp1 = bam_open(argv[1], "r");
    fp2 = bam_open(argv[2], "r");
    header1 = bam_header_read(fp1);
    header2 = bam_header_read(fp2);

    //Process the first file
    printf("reading in %s\n", argv[1]);
    while(bam_read1(fp1, read1) > 1) {
        if(++i + 100000 > max_length1) {
            max_length1 += 1000000;
            reads1 = realloc(reads1, sizeof(reads) * max_length1);
        }
        n1++;

        //id
        name = bam1_qname(read1);
        p = strchr(name, '_');
        if(p != NULL) *p = '\0';
//        p = strchr(name, '.');
//        (reads1+i)->id = atoi(++p);

        //tid
        (reads1+i)->tid = read1->core.tid;

        //pos1, pos2
        if(read1->core.flag & BAM_FPAIRED) {
            bam_read1(fp1, read2);
            if(read1->core.flag & BAM_FREAD1) {
                (reads1+i)->pos1 = read1->core.pos;
                (reads1+i)->pos2 = read2->core.pos;
            } else {
                (reads1+i)->pos2 = read1->core.pos;
                (reads1+i)->pos1 = read2->core.pos;
            }
        } else {
            (reads1+i)->pos1 = read1->core.pos;
            (reads1+i)->pos2 = 0;
        }
    }
    printf("\t%lu reads\n", i);

    //Process the second file
    printf("reading in %s\n", argv[2]);
    while(bam_read1(fp2, read1) > 1) {
        if(++j + 100000 > max_length2) {
            max_length2 += 1000000;
            reads2 = realloc(reads2, sizeof(reads) * max_length2);
        }
        n2++;

        //id
        name = bam1_qname(read1);
        p = strchr(name, '_');
        if(p != NULL) *p = '\0';
//        p = strchr(name, '.');
//        (reads2+j)->id = atoi(++p);

        //tid
        (reads2+j)->tid = read1->core.tid;

        //pos1, pos2
        if(read1->core.flag & BAM_FPAIRED) {
            bam_read1(fp2, read2);
            if(read1->core.flag & BAM_FREAD1) {
                (reads2+j)->pos1 = read1->core.pos;
                (reads2+j)->pos2 = read2->core.pos;
            } else {
                (reads2+j)->pos2 = read1->core.pos;
                (reads2+j)->pos1 = read2->core.pos;
            }
        } else {
            (reads2+j)->pos1 = read1->core.pos;
            (reads2+j)->pos2 = 0;
        }
    }
    printf("\t%lu reads\n", j);

    //Sort
    printf("Sorting %s\n", argv[1]);
    qsort(reads1, n1, sizeof(reads), cmp_reads);
    printf("Sorting %s\n", argv[2]);
    qsort(reads2, n2, sizeof(reads), cmp_reads);

    //Output the differences
    i=0;
    j=0;
    while(1) {
        if(i>=n1) break;
        if(j>=n2) break;

        rv = compare_reads(reads1+i, reads2+j);
        if(rv == 0) {
            i++;
            j++;
        } else if(rv == 1) {
            fprintf(of, "\t\t\t\t\t<\t%"PRId32"\t%s:%"PRId32"\t%s:%"PRId32"\n", \
                (reads2+j)->id, header2->target_name[(reads2+j)->tid], \
                (reads2+j)->pos1, header2->target_name[(reads2+j)->tid], (reads2+j)->pos2);
            j++;
        } else if(rv == -1) {
            fprintf(of, "%"PRId32"\t%s:%"PRId32"\t%s:%"PRId32"\t>\n", \
                (reads1+i)->id, header1->target_name[(reads1+i)->tid], \
                (reads1+i)->pos1, header1->target_name[(reads1+i)->tid], (reads1+i)->pos2);
            i++;
        } else {
            fprintf(of, "%"PRId32"\t%s:%"PRId32"\t%s:%"PRId32"\t!=", \
                (reads1+i)->id, header1->target_name[(reads1+i)->tid], \
                (reads1+i)->pos1, header1->target_name[(reads1+i)->tid], (reads1+i)->pos2);
            fprintf(of, "\t%"PRId32"\t%s:%"PRId32"\t%s:%"PRId32"\n", \
                (reads2+j)->id, header2->target_name[(reads2+j)->tid], \
                (reads2+j)->pos1, header2->target_name[(reads2+j)->tid], (reads2+j)->pos2);
            i++;
            j++;
        }
    }
    while(i<n1) {
        fprintf(of, "%"PRId32"\t%s:%"PRId32"\t%s:%"PRId32"\t>\n", \
            (reads1+i)->id, header1->target_name[(reads1+i)->tid], \
            (reads1+i)->pos1, header1->target_name[(reads1+i)->tid], (reads1+i)->pos2);
        i++;
    }
    while(j<n2) {
        fprintf(of, "\t\t\t\t\t<\t%"PRId32"\t%s:%"PRId32"\t%s:%"PRId32"\n", \
            (reads2+j)->id, header2->target_name[(reads2+j)->tid], \
            (reads2+j)->pos1, header2->target_name[(reads2+j)->tid], (reads2+j)->pos2);
        j++;
    }

    //Close up
    fclose(of);
    bam_close(fp1);
    bam_close(fp2);
    free(reads1);
    free(reads2);
    bam_header_destroy(header1);
    bam_header_destroy(header2);
    bam_destroy1(read1);
    bam_destroy1(read2);
    free(oname);
    return 0;
}
