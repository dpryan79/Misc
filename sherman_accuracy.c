#include "bam.h"

void usage(char *prog) {
    printf("Usage: %s file.bam\n", prog);
    printf("\n \
    This program will parse a Bison-aligned BAM file where reads were created\n \
    by Sherman (or any other read generator that tags reads with\n \
    X_chr?:start-stop coordinates as a name). Output is written to file.metrics\n \
    and consists of three columns: MAPQ, number correctly aligned, number\n \
    incorrectly aligned. Note that the BAM file MUST NOT be coordinate sorted\n \
    (Bison outputs a name sorted file)!\n \
\n \
-perfect  Ignore any alignments with a valid second best alignment (i.e. an XS\n \
    score).\n \
\n \
-all      Output by MAPQ instead of AS and XS score\n");
}

int main(int argc, char *argv[]) {
    char *oname, *p, *pstart, *location = malloc(1024*sizeof(char));
    char *XR, *XG;
    FILE *of = NULL;
    int start, stop, i, j, cstart, cstop, perfect = 0;
    int AS = 0, XS = 0;
    int correct[1000][1000], incorrect[1000][1000];
    int MAPQ_correct[47], MAPQ_incorrect[47];
    int all = 0;
    bamFile fp;
    bam_header_t *header;
    bam1_t *read1 = bam_init1();
    bam1_t *read2 = bam_init1();

    if(argc < 2) {
        usage(argv[0]);
        return 1;
    }
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            return 1;
        } else if(strcmp(argv[i], "-perfect") == 0) {
            perfect = 1;
        } else if(strcmp(argv[i], "-all") == 0) {
            all = 1;
        } else {
            fp = bam_open(argv[i], "r");
            oname = strdup(argv[i]);
            oname = realloc(oname, sizeof(char) * (strlen(oname) + 100));
            p = strrchr(oname, '.');
            *p = '\0';
            if(perfect) {
                sprintf(oname, "%s.no_secondary.metrics", oname);
            } else if(all) {
                sprintf(oname, "%s.all.metrics", oname);
            } else {
                sprintf(oname, "%s.metrics", oname);
            }
            of = fopen(oname, "w");
        }
    }

    //read in the header
    header = bam_header_read(fp);

    for(i=0; i<46; i++) {
        MAPQ_correct[i] = 0;
        MAPQ_incorrect[i] = 0;
    }
    for(i=0; i<1000; i++) {
        for(j=0; j<1000; j++) {
            correct[j][i] = 0;
            incorrect[j][i] = 0;
        }
    }

    //Loop over the reads
    while(bam_read1(fp, read1) > 1) {
        if(read1->core.flag & BAM_FPAIRED) {
            bam_read1(fp, read2);
            start = (read1->core.pos < read2->core.pos) ? read1->core.pos : read2->core.pos;
            stop = (bam_calend(&read1->core, bam1_cigar(read1)) > bam_calend(&read2->core, bam1_cigar(read2))) ? bam_calend(&read1->core, bam1_cigar(read1)) : bam_calend(&read2->core, bam1_cigar(read2));
        } else {
            start = read1->core.pos;
            stop = bam_calend(&read1->core, bam1_cigar(read1));
        }
        if(!all && (read1->core.qual == 0 || read1->core.flag & BAM_FUNMAP)) continue;

        //Do we care about reads with valid secondary alignments?
        XS = 0;
        if(perfect) {
            if(bam_aux_get(read1, "XS") != NULL) continue;
            AS = abs(bam_aux2i(bam_aux_get(read1, "AS")));
            if(read1->core.flag & BAM_FPAIRED) {
                if(bam_aux_get(read2, "XS") != NULL) continue;
                AS += abs(bam_aux2i(bam_aux_get(read2, "AS")));
            }
        } else {
            AS = abs(bam_aux2i(bam_aux_get(read1, "AS")));
            if(read1->core.flag & BAM_FPAIRED) {
                AS += abs(bam_aux2i(bam_aux_get(read2, "AS")));
            }
            if(bam_aux_get(read1, "XS") != NULL) {
                XS = abs(bam_aux2i(bam_aux_get(read1, "XS")));
            } else {
                if(read1->core.flag & BAM_FPAIRED) {
                    if(bam_aux_get(read2, "XS") != NULL) XS += abs(bam_aux2i(bam_aux_get(read1, "XS")));
                }
            }
        }

        //We need to adjust the start/stop coordinates depending on OT/OB/CTOT/CTOB strand
        XR = bam_aux2Z(bam_aux_get(read1, "XR"));
        XG = bam_aux2Z(bam_aux_get(read1, "XG"));
        if(strcmp(XR, "CT") == 0) {
            if(strcmp(XG, "CT") == 0) { //OT
                start++;
                stop++;
            }
        } else {
            if(strcmp(XG, "CT") == 0) { //CTOT
                start++;
                stop++;
            }
        }
        snprintf(location, 1024, "%s:%i-%i", header->target_name[read1->core.tid], start, stop);
        p = strchr(bam1_qname(read1), '_');
        p++;
        pstart = strtok(p, ":");
        cstart = atoi(strtok(NULL, "-"));
        cstop = atoi(strtok(NULL, "-"));
        //Allow a bit of wiggle room, since sometimes Sherman is off
        if(strncmp(p, header->target_name[read1->core.tid], strlen(header->target_name[read1->core.tid])) == 0 && abs(start-cstart) < 3 && abs(stop-cstop) < 3) {
            if(perfect) {
                correct[0][AS] += 1;
            } else if(all) {
                if(read1->core.qual > 45) read1->core.qual=45;
                MAPQ_correct[read1->core.qual] += 1;
            } else {
                correct[XS][AS] += 1;
            }
        } else {
            if(perfect) {
                incorrect[0][AS] += 1;
            } else if(all) {
                if(read1->core.qual > 45) read1->core.qual=45;
                MAPQ_incorrect[read1->core.qual] += 1;
            } else {
                incorrect[XS][AS] += 1;
            }
        }
    }

    if(perfect) {
        fprintf(of, "AS\tcorrect\tincorrect\n");
        for(i=0; i<1000; i++) {
            if(correct[0][i] + incorrect[0][i] > 0) {
                fprintf(of, "%i\t%i\t%i\n", -1*i, correct[0][i], incorrect[0][i]);
            }
        }
    } else if(all) {
        fprintf(of, "MAPQ\tcorrect\tincorrect\n");
        for(i=0; i<46; i++) {
            fprintf(of, "%i\t%i\t%i\n", i, MAPQ_correct[i], MAPQ_incorrect[i]);
        }
    } else {
        fprintf(of, "AS\tXS\tcorrect\tincorrect\n");
        for(i=0; i<1000; i++) {
            for(j=0; j<1000; j++) {
                if(correct[j][i] + incorrect[j][i] > 0) {
                    fprintf(of, "%i\t%i\t%i\t%i\n", -1*i,-1*j, correct[j][i], incorrect[j][i]);
                }
            }
        }
    }

    //Clean up
    bam_destroy1(read1);
    bam_destroy1(read2);
    bam_header_destroy(header);
    bam_close(fp);
    fclose(of);

    return 0;
}
