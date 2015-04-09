#include "../htslib/htslib/sam.h"
#include <inttypes.h>
#include <pthread.h>

typedef struct {
    int32_t tid;
    int32_t start;
    int32_t end;
} QueuedRegion;

struct QueuedRegions {
    int32_t nJobs;
    QueuedRegion **regions;
} QueuedRegions;

pthread_mutex_t mutQueueJobNum;
int32_t QueueJobNum;
struct QueuedRegions Regions;

char *fileName;
uint32_t *counts;

void * WorkerThread(void *a);


void queue_init();
void queue_destroy();
void queue_push(QueuedRegion *reg);


