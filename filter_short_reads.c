//gcc -Wall -O3 -o ~/bin/filter_short_reads -I/home/ryand/include -L/home/ryand/lib filter_short_reads.c -lhts -lz -lpthread -lm
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>

//take a stack of multimappers originating from the same read and determine if they can all be assigned to a single transcript
void process_stack(bam1_t **reads, int nreads, bam_hdr_t *header, int max_read_length, htsFile *ofile) {
    int i = 0, best_hit = -1;
    int32_t start, end;
    int *good_match, *possible_match;
    uint32_t cigar;
    uint8_t nfound = 0, j, possible_found = 0;
    double diff; //The percentage difference in the amount of a transcript covered and the transcript's length

    if(reads[0]->core.l_qseq > 0.9*max_read_length) { //The read is too long to accurately assign
        possible_match = (int*) calloc(nreads, sizeof(int));
        for(i=0;i<nreads;i++) {
            if((reads[i]->core.flag) & 16) continue; //Ignore reversed reads
            possible_match[i] = 1; //No reverse complemented, so it's possibly correct
            possible_found++;
        }
        best_hit = -1;
        for(i=0, j=1; i<nreads; i++) {
            if(!(possible_match[i])) continue;
            if(best_hit == -1) {
                best_hit = i;
                if(reads[i]->core.flag & 256) reads[i]->core.flag -= 256;
            } else {
                reads[i]->core.flag |= 256;
            }
            //replace the current NH and HI tags
            bam_aux_del(reads[i], bam_aux_get(reads[i], "NH"));
            bam_aux_del(reads[i], bam_aux_get(reads[i], "HI"));
            bam_aux_append(reads[i], "NH", 'c', 1, &nfound); //Replace NH
            bam_aux_append(reads[i], "HI", 'c', 1, &j); //Replace HI
            bam_write1(ofile->fp.bgzf, reads[i]);
            j++;
        }
        for(i=0;i<nreads;i++) bam_destroy1(reads[i]);
        free(possible_match);
    } else {
        good_match = (int*) calloc(nreads, sizeof(int));
        possible_match = (int*) calloc(nreads, sizeof(int));
        for(i=0;i<nreads;i++) {
            if((reads[i]->core.flag) & 16) continue; //Ignore reversed reads
            possible_match[i] = 1; //No reverse complemented, so it's possibly correct
            possible_found++;
            start = reads[i]->core.pos;
            cigar = *bam_get_cigar(reads[i]);
            if((cigar & 0xF) == BAM_CSOFT_CLIP) start -= (cigar >> 4); //Deal with soft-clipping
            end = bam_endpos(reads[i]);
            //Deal with 3' soft clipping
            cigar = *(bam_get_cigar(reads[i]) + reads[i]->core.n_cigar - 1);
            if((cigar & 0xF) == BAM_CSOFT_CLIP) end += (cigar >> 4);
            diff = fabs(1.0 - (end-start)/((double)header->target_len[reads[i]->core.tid]));
            if(diff <= 0.05) good_match[i] = 1;
        }
        //Do we have one target with perfect occupancy?
        for(i=0; i<nreads; i++) {
            if(good_match[i]) {
                if(best_hit<0) best_hit = i; //Keep the best_hit as the first good hit
                nfound++;
            }
        }
        if(nfound) {
            for(i=best_hit, j=1; i<nreads; i++) {
                if(!good_match[i]) continue;
                if(i==best_hit) {
                    if(reads[i]->core.flag & 256) reads[i]->core.flag -= 256;
                } else {
                    reads[i]->core.flag |= 256;
                }
                //replace the current NH and HI tags
                bam_aux_del(reads[i], bam_aux_get(reads[i], "NH"));
                bam_aux_del(reads[i], bam_aux_get(reads[i], "HI"));
                bam_aux_append(reads[i], "NH", 'c', 1, &nfound); //Replace NH
                bam_aux_append(reads[i], "HI", 'c', 1, &j); //Replace HI
                bam_write1(ofile->fp.bgzf, reads[i]);
                j++;
            }
        } else {
            best_hit = -1;
            for(i=0, j=1; i<nreads; i++) {
                if(!(possible_match[i])) continue;
                if(best_hit == -1) {
                    best_hit = i;
                    if(reads[i]->core.flag & 256) reads[i]->core.flag -= 256;
                } else {
                    reads[i]->core.flag |= 256;
                }
                //replace the current NH and HI tags
                bam_aux_del(reads[i], bam_aux_get(reads[i], "NH"));
                bam_aux_del(reads[i], bam_aux_get(reads[i], "HI"));
                bam_aux_append(reads[i], "NH", 'c', 1, &possible_found); //Replace NH
                bam_aux_append(reads[i], "HI", 'c', 1, &j); //Replace HI
                bam_write1(ofile->fp.bgzf, reads[i]);
                j++;
            }
        }

        for(i=0; i<nreads; i++) bam_destroy1(reads[i]);
        free(good_match);
        free(possible_match);
    }
}

void usage(char *prog) {
    printf("Usage: samtools view -h file.bam | %s [OPTIONS] output.bam\n", prog);
    printf("\n\
The problem:\n\
Reads originating from short RNAs often map equally well to multiple\n\
transcripts, but we would ideally like to assign them to only one transcript. To\n\
do so, we can utilize the fact that reads are often longer than the transcripts\n\
(even after trimming).\n\
\n\
An example:\n\
Suppose a read aligns equally well to three transcripts.\n\
\n\
|+++++++++++++++++++++++++++++???????????????| The original read\n\
|+++++++++++++++++++++++++++++| The trimmed read\n\
|--------------------------------| Transcript A\n\
|-------------------------------------| Transcript B\n\
|-----------------------------| Transcript C\n\
\n\
Given how trimming works, it's most likely tha the read arose from transcript C.\n\
We will asign a read to a single transcript if and only if its bounds exactly\n\
match that of a single transcript to which it can match. Further, the reads\n\
length can't be more than 90%% of it's untrimmed length. This is to limit false-\n\
positive assignments due to a ambiguity about whether the read actually contains\n\
the entirety of the sequence of the transcript from which it arose. Thus, the\n\
following example would not be assigned to a transcript:\n\
\n\
|+++++++++++++++++++++++++++++?| The original read\n\
|+++++++++++++++++++++++++++++| The trimmed read\n\
|--------------------------------| Transcript A\n\
|-------------------------------------| Transcript B\n\
|-----------------------------| Transcript C\n\
\n\
It should be noted that an assignment will not be made to a transcript it if\n\
requires reverse-complement mapping, due to the strand-specific nature of the\n\
library-prep.\n\
\n\
Options\n\
\n\
-L	The original read length (default 50).\n\
\n\
-@	Number of compression threads (default 1).\n\
\n\
--keep-multi	If set, then reads that can't be assigned to a single transcript\n\
	will be output as is rather than being discarded.\n");
}

int main(int argc, char *argv[]) {
    bam_hdr_t *header = NULL;
    bam1_t **reads = NULL, *read = NULL;
    htsFile *ifile, *ofile = NULL;
    int length = 50, i, m_stack, NH, nthreads = 1;
    uint8_t *p;
    char *oname = NULL;
    int rv = 0;

    if(argc==1) {
        usage(argv[0]);
        return 0;
    }
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else if(strcmp(argv[i], "-L") == 0) {
            length = atoi(argv[++i]);
        } else if(strcmp(argv[i], "-@") == 0) {
            nthreads = atoi(argv[++i]);
        } else if(oname == NULL) {
            oname = argv[i];
        } else {
            usage(argv[0]);
            return -1;
        }
    }

    if(oname == NULL) {
        usage(argv[0]);
        return -1;
    }
    ofile = sam_open(oname, "wb");

    if(nthreads > 1) bgzf_mt(ofile->fp.bgzf, nthreads, 256);

    ifile = sam_open("-", "r");
    header = sam_hdr_read(ifile);
    bam_hdr_write(ofile->fp.bgzf, header);

    m_stack = 100;
    read = bam_init1();
    reads = malloc(sizeof(bam1_t*) * m_stack);
    while(sam_read1(ifile, header, read) >= 0) {
        p = bam_aux_get(read, "NH");
        NH = 1;
        if(p != NULL) NH = bam_aux2i(p);
        if(NH == 1 && (read->core.flag & 16) == 0) {
            bam_write1(ofile->fp.bgzf, read);
        } else {
            if(NH > m_stack) {
                reads = realloc(reads, sizeof(bam1_t *) * NH);
                if(reads == NULL) {
                    printf("Couldn't reallocate enough room to hold %i reads. Perhaps there's not enough memory!", NH);
                    rv = -1;
                    goto err;
                }
            }
            *reads = bam_dup1(read);
            for(i=1; i<NH; i++) {
                if(sam_read1(ifile, header, read) < 0) {
                    printf("There was an error reading the input, or we've reached the end of the file earlier than expected.");
                    rv = -1;
                    goto err;
                }
                *(reads+i) = bam_dup1(read);
            }
            process_stack(reads, NH, header, length, ofile);
        }
    }

    free(reads);
err:
    bam_destroy1(read);
    bam_hdr_destroy(header);
    sam_close(ofile);
    sam_close(ifile);

    return rv;
}
