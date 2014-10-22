//gcc -Wall -o ~/bin/shortReads2BAMandCounts -I/home/ryand/include shortReads2BAMandCounts.c /home/ryand/lib/libhts.a -lz -lpthread
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

int32_t parse_mappings(FILE *fp, char ***transcript_names, int32_t **tsid2cid, int32_t *ntranscripts) {
    char *line = malloc(sizeof(char)*1024); //This should suffice
    char *p;
    int32_t nclusters = 0;
    *ntranscripts = 0;

    while(fgets(line, 1024, fp)) {
        p = strtok(line, "\t");
        *transcript_names = realloc(*transcript_names, (1+*ntranscripts) * sizeof(char*));
        *tsid2cid = realloc(*tsid2cid, (1+*ntranscripts) * sizeof(int32_t));
        (*transcript_names)[*ntranscripts] = strdup(p); //This must be free()d
        p = strtok(NULL, "\n");
        (*tsid2cid)[*ntranscripts] = strtoll(p, NULL, 10);
        if(nclusters < strtoll(p, NULL, 10)) nclusters = strtoll(p, NULL, 10);
        *ntranscripts += 1;
    }
    free(line);

    return nclusters;
}

//Produce a map of hdr->target_name[i] to transcript_names[j]
int32_t *mapBAM2Cluster(bam_hdr_t *hdr, char **transcript_names, int32_t *tsid2cid) {
    int32_t i, j;
    int32_t *BAM2Cluster = malloc(sizeof(int32_t) * hdr->n_targets);
    if(BAM2Cluster == NULL) {
        fprintf(stderr, "Couldn't allocate enough space for BAM2Cluster!!!\n");
        return NULL;
    }

    for(i=0; i<hdr->n_targets; i++) {
        if(strcmp(hdr->target_name[i], transcript_names[i]) == 0) BAM2Cluster[i] = tsid2cid[i];
        else {
            BAM2Cluster[i] = -1;
            for(j=0; j<hdr->n_targets; j++) {
                if(strcmp(hdr->target_name[i], transcript_names[j]) == 0) {
                    BAM2Cluster[i] = tsid2cid[j];
                    break;
                }
            }
            if(BAM2Cluster[i] == -1) {
                fprintf(stderr, "Couldn't find a match for %s in the RNA2Cluster file!!!\n", hdr->target_name[i]);
                return NULL;
            }
        }
    }

    return BAM2Cluster;
}

//take a stack of multimappers originating from the same read and determine if they can all be assigned to a single transcript
//return -1 on error!
int process_stack(bam1_t **reads, int nreads, bam_hdr_t *header, int32_t *BAM2Cluster, int max_read_length, htsFile *ofile, uint32_t *counts) {
    int i = 0, best_hit = -1, nfound = 0, rv = 0, same_cluster = 1, npossible = 0;
    int32_t start, end;
    uint32_t *diff, *possible; //This should be overkill
    uint32_t cigar;

    possible = (uint32_t*) calloc(nreads, sizeof(uint32_t)); //Need to check for NULL!!
    if(reads[0]->core.l_qseq > max_read_length-10) { //The read is too long to accurately restrict
        for(i=0;i<nreads;i++) {
            if(reads[i]->core.flag & 0xf) continue; //Ignore reversed reads
            if(best_hit < 0) best_hit = i;
            possible[i] = 1;
            nfound++;
        }
        if(nfound) {
            for(i=0; i<nreads; i++) {
                if(!possible[i]) continue;
                rv++;
                //replace the current NH and HI tags
                bam_aux_del(reads[i], bam_aux_get(reads[i], "NH"));
                bam_aux_del(reads[i], bam_aux_get(reads[i], "HI"));
                bam_aux_append(reads[i], "NH", 'c', 1, (uint8_t *) &nfound); //Replace NH
                bam_aux_append(reads[i], "HI", 'c', 1, (uint8_t *) &rv); //Replace HI
                //Append an XC:i tag, denoting the cluster membership
                bam_aux_append(reads[i], "XC", 'i', 4, (uint8_t*) &(BAM2Cluster[reads[i]->core.tid])); //Append XC
                //fix the flag
                if(rv==1) {if(reads[i]->core.flag & 0x100) reads[i]->core.flag -= 0x100; }
                else reads[i]->core.flag |= 0x100;
                bam_write1(ofile->fp.bgzf, reads[i]);

                if(BAM2Cluster[reads[best_hit]->core.tid] != BAM2Cluster[reads[i]->core.tid]) same_cluster = 0;
            }
            rv = 0;
            if(same_cluster) counts[BAM2Cluster[reads[best_hit]->core.tid]] += 1;
        }
    } else {
        diff = (uint32_t*) calloc(nreads, sizeof(uint32_t)); //Need to check for NULL!!
        for(i=0;i<nreads;i++) {
            if(reads[i]->core.flag & 0xf) continue; //Ignore reversed reads
            possible[i] = 1;
            npossible++;
            start = reads[i]->core.pos;
            cigar = *bam_get_cigar(reads[i]);
            if((cigar & 0xF) == BAM_CSOFT_CLIP) start -= (cigar >> 4); //Deal with soft-clipping
            end = bam_endpos(reads[i]);
            //Deal with 3' soft clipping
            cigar = *(bam_get_cigar(reads[i]) + reads[i]->core.n_cigar - 1);
            if((cigar & 0xF) == BAM_CSOFT_CLIP) end += (cigar >> 4);
            diff[i] = abs(end-start-header->target_len[reads[i]->core.tid]);
        }
        /***********************************************************************
        * Find the alignment with the smallest difference between mapped and
        * target length. Only accept differences <= 6. If there are multiple
        * equally best hits, then we'll first see if they belong to the same
        * cluster.
        ***********************************************************************/
        for(i=0; i<nreads; i++) {
            if(diff[i] < 7) {
                if(best_hit == -1) {
                    best_hit = i;
                    nfound++;
                } else {
                    if(diff[i] < diff[best_hit]) {
                        best_hit = i;
                        nfound = 1;
                    } else if(diff[i] == diff[best_hit]) nfound++;
                }
            }
        }
        if(nfound) {
            for(i=best_hit; i<nreads; i++) {
                if(diff[i] > diff[best_hit]) continue; //Only output the top stratum of hits
                if(BAM2Cluster[reads[best_hit]->core.tid] != BAM2Cluster[reads[i]->core.tid]) same_cluster = 0;
                rv++;
                //replace the current NH and HI tags
                bam_aux_del(reads[i], bam_aux_get(reads[i], "NH"));
                bam_aux_del(reads[i], bam_aux_get(reads[i], "HI"));
                bam_aux_append(reads[i], "NH", 'c', 1, (uint8_t *) &nfound); //Replace NH
                bam_aux_append(reads[i], "HI", 'c', 1, (uint8_t *) &rv); //Replace HI
                //Append an XC:i tag, denoting the cluster membership
                bam_aux_append(reads[i], "XC", 'i', 4, (uint8_t*) &(BAM2Cluster[reads[i]->core.tid])); //Append XC
                //fix the flag
                if(i == best_hit) { if(reads[i]->core.flag & 0x100) reads[i]->core.flag -= 0x100; }
                else reads[i]->core.flag |= 0x100;
                bam_write1(ofile->fp.bgzf, reads[i]);
            }
            if(same_cluster) counts[BAM2Cluster[reads[best_hit]->core.tid]] += 1;
        } else { //No acceptable hits
            for(i=0; i<nreads; i++) if(possible[i]) nfound++; //repurposing the variable
            for(i=0; i<nreads; i++) {
                if(possible[i] == 0) continue;
                if(!rv) { if(reads[i]->core.flag & 0x100) reads[i]->core.flag -= 0x100; }
                else reads[i]->core.flag |= 0x100;
                rv++;
                //replace the current NH and HI tags
                bam_aux_del(reads[i], bam_aux_get(reads[i], "NH"));
                bam_aux_del(reads[i], bam_aux_get(reads[i], "HI"));
                bam_aux_append(reads[i], "NH", 'c', 1, (uint8_t*) &nfound); //Replace NH
                bam_aux_append(reads[i], "HI", 'c', 1, (uint8_t*) &rv); //Replace HI
                //Append an XC:i tag, denoting the cluster membership
                bam_aux_append(reads[i], "XC", 'i', 4, (uint8_t *) &(BAM2Cluster[reads[i]->core.tid])); //Append XC
                bam_write1(ofile->fp.bgzf, reads[i]);
            }
            rv = 0; //don't mislead the rest of the program!
        }
        free(diff);
    }
    free(possible);
    return rv;
}

void usage(char *prog) {
    printf("Usage: samtools view -h file.bam | %s [OPTIONS] RNA2Cluster.txt output_prefix\n", prog);
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
-@	Number of compression threads (default 1).\n");
}

int main(int argc, char *argv[]) {
    bam_hdr_t *header = NULL;
    bam1_t **reads = NULL, *read = NULL;
    htsFile *ifile, *obam= NULL;
    FILE *otxt, *mapping_file;
    int length = 50, i, m_stack = 100, NH, rv = 0, nthreads = 1;
    uint8_t *p;
    uint32_t original_count = 0, new_count = 0;
    char *oname_prefix = NULL, *oname, *mapping_file_name = NULL;
    char **transcript_names = NULL;
    uint32_t *cluster_counts = NULL; //Overkill? Why yes of course!
    int32_t *BAM2Cluster, nclusters = 0, ntranscripts = 0;
    int32_t *tsid2cid = NULL;

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
        } else if(mapping_file_name == NULL) {
            mapping_file_name = argv[i];
        } else if(oname_prefix == NULL) {
            oname_prefix = argv[i];
        } else {
            usage(argv[0]);
            return -1;
        }
    }

    if(oname_prefix == NULL) {
        usage(argv[0]);
        return -1;
    }

    printf("Writing alignments to %s.bam and counts to %s.txt\n", oname_prefix, oname_prefix);
    oname = malloc(sizeof(char) * (strlen(oname_prefix) + 5));
    sprintf(oname, "%s.bam", oname_prefix);
    obam = sam_open(oname, "wb");
    sprintf(oname, "%s.txt", oname_prefix);
    otxt = fopen(oname, "w");
    free(oname);
    mapping_file = fopen(mapping_file_name, "r");

    //Read in the transcript->cluster mapping file
    nclusters = parse_mappings(mapping_file, &transcript_names, &tsid2cid, &ntranscripts); //ntranscripts is set by this too
    fclose(mapping_file);

    if(nthreads > 1) bgzf_mt(obam->fp.bgzf, nthreads, 256);

    //Open the input
    ifile = sam_open("-", "r");
    header = sam_hdr_read(ifile);
    bam_hdr_write(obam->fp.bgzf, header);
    if(ntranscripts != header->n_targets) {
        fprintf(stderr, "There's a mismatch in the number of transcripts between the input SAM file (%" PRId32 ") and %s (%" PRId32 ").\n", header->n_targets, mapping_file_name, ntranscripts);
        rv = -3;
        goto err;
    }
    cluster_counts = calloc(nclusters, sizeof(uint32_t));

    //Map the header transcript names to the mapping file transcript names
    BAM2Cluster = mapBAM2Cluster(header, transcript_names, tsid2cid);

    //Allocate the stack that'll hold the mappings
    m_stack = 100;
    reads = malloc(sizeof(bam1_t*) * m_stack);
    if(reads == NULL) {
        fprintf(stderr, "Couldn't allocate enough room for the 'reads' array. Out of memory?\n");
        rv =-2;
        goto err;
    }
    for(i=0; i<m_stack; i++) reads[i] = bam_init1();

    read = bam_init1();
    while(sam_read1(ifile, header, read) >= 0) {
        p = bam_aux_get(read, "NH");
        NH = 1;
        if(p != NULL) NH = bam_aux2i(p);
        original_count++;
        //Grow the stack as needed
        if(NH > m_stack) {
            if(realloc(reads, sizeof(bam1_t *) * NH) == NULL) {
                fprintf(stderr, "Couldn't reallocate enough room to hold %i reads. Perhaps there's not enough memory!", NH);
                rv = -1;
                goto err;
            }
            for(i=m_stack; i<NH; i++) reads[i] = bam_init1();
            m_stack = NH;
        }
        reads[0] = bam_copy1(reads[0], read);
        for(i=1; i<NH; i++) {
            if(sam_read1(ifile, header, read) < 0) {
                fprintf(stderr, "There was an error reading the input, or we've reached the end of the file earlier than expected.");
                rv = -1;
                goto err;
            }
            reads[i] = bam_copy1(reads[i], read);
        }
        new_count += process_stack(reads, NH, header, BAM2Cluster, length, obam, cluster_counts);
    }

    fprintf(otxt, "ClusterID\t%s\n", oname_prefix);
    for(i=0; i<nclusters; i++) fprintf(otxt, "cluster %i\t%" PRIu32 "\n", i+1, cluster_counts[i]);

err:
    for(i=0; i<header->n_targets; i++) free(transcript_names[i]);
    free(transcript_names);
    for(i=0; i<m_stack; i++) bam_destroy1(reads[i]);
    free(reads);
    if(cluster_counts) free(cluster_counts);
    fclose(otxt);
    bam_destroy1(read);
    sam_close(obam);
    sam_close(ifile);
    free(tsid2cid);

    printf("Would have originally reported %" PRIu32 " alignments and have now reported %" PRIu32 "\n", original_count, new_count);
    return rv;
}
