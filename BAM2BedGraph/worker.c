#include <pthread.h>
#include <stdio.h>
#include "../htslib/htslib/sam.h"
#include "../htslib/htslib/hts.h"
#include "BAM2BedGraph.h"

//Determine the bounds of a fragment, ignoring 
//This isn't ideal, we ignore for PE datasets, we ignore singletons
int getbounds(bam1_t *b, bam_hdr_t *hdr, int32_t *left, int32_t *right) {
    if(b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FDUP)) return 0;
    if(b->core.flag & BAM_FREAD2) return 0;
    if(b->core.flag & BAM_FPAIRED) {
        if(!(b->core.flag & BAM_FPROPER_PAIR)) return 0;
    }

    if(b->core.flag & 16) {
        *right = bam_endpos(b);
        if(b->core.flag&2) {
            *left = *right - b->core.isize;
        } else {
            *left = *right - 300;
        }
    } else {
        *left = b->core.pos;
        if(b->core.flag & 2) {
            *right = b->core.pos + b->core.isize;
         } else {
             *right = b->core.pos + 300;
         }
    }
    if(*left < 0) *left = 0;
    if(*right >= hdr->target_len[b->core.tid]) *right = hdr->target_len[b->core.tid];
    return 1;
}

/*
 * Returns a NULL pointer, regardless of error or not (this is just a test program)
 * The input is actually ignored.
 */
void * WorkerThread(void *a) {
    QueuedRegion *region;
    int32_t jobID, lbound, rbound, i;
    htsFile *fp;
    hts_idx_t *idx;
    hts_itr_t *iter;
    bam_hdr_t *hdr;
    bam1_t *b;

    while(1) {
        //get the region
        pthread_mutex_lock(&mutQueueJobNum);
        jobID = QueueJobNum++;
        pthread_mutex_unlock(&mutQueueJobNum);
        if(jobID >= Regions.nJobs) return NULL;
        region = Regions.regions[jobID];
fprintf(stderr, "Thread %"PRId32" got region %"PRId32" %"PRId32"-%"PRId32"\n", *((int32_t*) a), region->tid, region->start, region->end);

        fp = sam_open(fileName, "r");
        idx = sam_index_load(fp, fileName);
        hdr = sam_hdr_read(fp);
        b = bam_init1();

        //Iterate
        iter = sam_itr_queryi(idx, region->tid, region->start, region->end);
        while(sam_itr_next(fp, iter, b) >= 0) {
            if(getbounds(b, hdr, &lbound, &rbound)) {
                for(i=lbound; i<=rbound; i++) {
                    counts[i]++;
                }
            }
        }
        sam_itr_destroy(iter);
        bam_destroy1(b);
        sam_close(fp);
        bam_hdr_destroy(hdr);
        hts_idx_destroy(idx);
    }
    return NULL;
}
