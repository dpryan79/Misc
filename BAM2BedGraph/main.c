#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "../htslib/htslib/sam.h"
#include "BAM2BedGraph.h"

void dumpCounts(char *chrom, int32_t l) {
    int32_t i;
    for(i=0; i<l; i++) {
        if(!counts[i]) continue;
        printf("%s\t%"PRId32"\t%"PRId32"\t%"PRIu32"\n", chrom, i, i+1, counts[i]);
    }
}

void usage() {
    printf("Usage: BAM2BedGraph <input>\n");
    printf("\n"
"This tools accepts a local or remote indexed BAM/CRAM file and produces a\n"
"BedGraph file with approximate per-base depth. This is meant as a test program.\n"
"4 worker threads will work on this, a la deeptools. Output is to the console.\n"
);
}

int main(int argc, char *argv[]) {
    htsFile *fp = NULL;
    bam_hdr_t *hdr = NULL;
    int32_t i, j, IDs[4] = {0, 1, 2, 3};
    int rv;
    QueuedRegion *reg;
    pthread_t threads[4];
    
    if(argc != 2 || strcmp("-h", argv[1]) == 0) {
        usage();
        return 0;
    }

    fp = sam_open(argv[1], "r");
    if(!fp) {
        fprintf(stderr, "Couldn't open %s!\n", argv[1]);
        return 1;
    }
    hdr = sam_hdr_read(fp);
    if(!hdr) {
        fprintf(stderr, "Couldn't read in the header!\n");
        sam_close(fp);
        return 1;
    }

    fileName = argv[1];
    for(i=0; i<hdr->n_targets; i++) {
        //Initialize the globals
        counts = calloc(hdr->target_len[i], sizeof(uint32_t));
        assert(counts);
        queue_init();

        j = 0;
        for(j=0; j < hdr->target_len[i]; j+=1e5) {
            reg = malloc(sizeof(QueuedRegion));
            assert(reg);
            reg->tid = i;
            reg->start = j;
            reg->end = 1e5+j+1; //Is this 1 or 0-based?
            queue_push(reg);
        }

        //start the workers
        for(j=0; j<4; j++) {
            rv = pthread_create(&(threads[j]), NULL, WorkerThread, (void*) &IDs[j]);
            if(rv) fprintf(stderr, "pthread_create returned %i\n", rv);
        }
        for(j=0; j<4; j++) pthread_join(threads[j], NULL);

        //print the output
        dumpCounts(hdr->target_name[i], hdr->target_len[i]);

        //clean up
        queue_destroy();
        free(counts);
    }

    bam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}
