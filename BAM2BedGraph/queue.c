#include "BAM2BedGraph.h"
#include <stdlib.h>
#include <assert.h>

void queue_init() {
    Regions.nJobs = 0;
    Regions.regions = NULL;
}

void queue_destroy() {
    int32_t i;
    if(Regions.nJobs) {
        for(i=0; i<Regions.nJobs; i++) {
            free(Regions.regions[i]);
        }
    }
    if(Regions.regions) free(Regions.regions);
    Regions.nJobs = 0;
    QueueJobNum = 0;
}

//This is terribly inefficient and should never be used in a real program!
void queue_push(QueuedRegion *reg) {
    Regions.nJobs++;
    Regions.regions = realloc(Regions.regions, Regions.nJobs*sizeof(QueuedRegion*));
    assert(Regions.regions);
    Regions.regions[Regions.nJobs-1] = reg;
}
