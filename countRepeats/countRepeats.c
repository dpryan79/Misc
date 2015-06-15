//gcc -g -L/home/ryan/Documents/SourceCode/libGTF -I/home/ryan/Documents/SourceCode/libGTF -L/home/ryan/Documents/SourceCode/htslib -I/home/ryan/Documents/SourceCode/htslib -o getChunk getChunk.c -lhts -lGTF -lz -lpthread -lpcre
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <pthread.h>
#include "htslib/sam.h"
#include "gtf.h"

uint64_t uniqueName, dupName;
uint64_t uniqueClass, dupClass;
uint64_t uniqueFamily, dupFamily;
cntTable *g_cntName, *g_cntFamily, *g_cntClass;

int exons(void *l) {
    GTFline *line = (GTFline *) l;
    if(strcmp(line->feature.s, "exon") == 0) return 1;
    return 0;
}

int passThrough(void *l) {
    return 1;
}

int passThroughEntry(GTFtree *t, GTFentry *e) {
    return 1;
}

int cmpBAMentries(GTFentry *e0, GTFentry *e1) {
    return 1;
}

typedef struct {
    int l,m;
    bam1_t **b;
} aList;

aList *al_init() {
    aList *al = calloc(1, sizeof(aList));
    assert(al);
    return al;
}

void al_reset(aList *al) {
    int i;
    for(i=0; i<al->l; i++) bam_destroy1(al->b[i]);
    al->l = 0;
}

void al_destroy(aList *al) {
    al_reset(al);
    if(al->m) free(al->b);
    free(al);
}

void al_push(aList *al, bam1_t *b) {
    int i;
    if(al->l+1 >= al->m) {
        al->m += 32;
        al->b = realloc(al->b, sizeof(bam1_t*)*(al->m));
        assert(al->b);
        for(i=al->l; i<al->m; i++) al->b[i] = NULL;
    }
    al->b[al->l] = bam_dup1(b);
    al->l++;
}

int al_sort_func(const void *a, const void *b) {
    bam1_t *pa = *(bam1_t**)a;
    bam1_t *pb = *(bam1_t**)b;
    int AS1, AS2;
    uint8_t *auxp;
    auxp = bam_aux_get(pa, "AS");
    AS1 = bam_aux2i(auxp);
    auxp = bam_aux_get(pb, "AS");
    AS2 = bam_aux2i(auxp);
    if(AS1 > AS2) return -1;
    else if(AS1 < AS2) return 1;
    return 0;
}

void al_sort(aList *al) {
    qsort((void *) al->b, al->l, sizeof(bam1_t**), al_sort_func);
}
        
aList *getGroup(aList *al, samFile *fp, bam_hdr_t *hdr) {
    bam1_t *b = bam_init1();
    int NH = 1, i;
    uint8_t *auxp;
    al_reset(al);

    if(sam_read1(fp, hdr, b) < 0) {
        al_destroy(al);
        bam_destroy1(b);
        return NULL;
    }
    auxp = bam_aux_get(b, "NH");
    if(auxp) NH = bam_aux2i(auxp);
    al_push(al, b);
    for(i=1; i<NH; i++) {
        assert(sam_read1(fp, hdr, b) >= 0);
        al_push(al, b);
    }
    bam_destroy1(b);
    return al;
}

int cntTop(aList *al) {
    int i, n = 1, AS;
    AS = bam_aux2i(bam_aux_get(al->b[0], "AS"));
    for(i=1; i<al->l; i++) {
        if(AS > bam_aux2i(bam_aux_get(al->b[i], "AS"))) break;
        n++;
    }
    return n;
}

//getGroup
//Possibly al_sort the cntTop
//The overlapsGenes, possibly al_reset, then overlapsRepeats
//Print a repeat family/etc. accordingly
int overlapsGenes(GTFtree *t, bam_hdr_t *hdr, aList *al, int l, int strandType, int matchType) {
    int i;
    overlapSet *os;
    for(i=0; i<l; i++) {
        os = findOverlapsBAM(t, al->b[i], hdr, matchType, strandType, NULL, cmpBAMentries);
        if(os->l) break;
    }
    if(os->l) {
        os_destroy(os);
        return 1;
    } else {
        os_destroy(os);
        return 0;
    }
}

void overlapsRepeats(GTFtree *t, bam_hdr_t *hdr, aList *al, int l, int strandType, int matchType, cntTable *cntName, cntTable *cntClass, cntTable *cntFamily) {
    int i;
    overlapSet *os;
    overlapSetList *osl = osl_init();
    uniqueSet *family, *class, *name;
    for(i=0; i<l; i++) {
        os = findOverlapsBAM(t, al->b[i], hdr, matchType, strandType, passThroughEntry, cmpBAMentries);
        if(os->l) osl_push(osl, os);
        else os_destroy(os);
    }
    //Take the union
    os = osl_union(osl);
    name = uniqueAttributes(os, "repName");
    class = uniqueAttributes(os, "repClass");
    family = uniqueAttributes(os, "repFamily");

    //Update the counts (should do some expectation maximization...)
    if(name) {
        if(name->l == 1) {
            incCntTable(cntName, us_val(name, 0));
        }
    }
    if(class) {
        if(class->l == 1) {
            incCntTable(cntClass, us_val(class, 0));
        }
    }
    if(family) {
        if(family->l == 1) {
            incCntTable(cntFamily, us_val(family, 0));
        }
    }

    //Clean up
    if(os) os_destroy(os);
    if(osl) osl_destroy(osl);
    if(name) us_destroy(name);
    if(class) us_destroy(class);
    if(family) us_destroy(family);
}

typedef struct {
    GTFtree *genes;
    GTFtree *repeats;
    htsFile *fp;
    bam_hdr_t *hdr;
} worker_in;

pthread_mutex_t fp_lock;
pthread_mutex_t num_lock;

void mergeCntTable(cntTable *g, cntTable *a) {
    int32_t i;
    for(i=0; i<g->ht->l; i++) g->cnts[i] += a->cnts[i];
}

void printCntTable(cntTable *ct, char *type) {
    int32_t i;
    for(i=0; i<ct->ht->l; i++) printf("%s\t%"PRIu32"\t%s\n", ct->ht->str[i], ct->cnts[i], type);
}

static void *worker(void *p) {
    worker_in *win = (worker_in*) p;
    aList *al = al_init();
    uint64_t l_uniqueName=0, l_dupName=0;
    uint64_t l_uniqueClass=0, l_dupClass=0;
    uint64_t l_uniqueFamily=0, l_dupFamily=0;
    cntTable *cntName, *cntFamily, *cntClass;
    GTFtree *genes = win->genes;
    GTFtree *repeats = win->repeats;
    htsFile *fp = win->fp;
    bam_hdr_t *hdr = win->hdr;
    cntName = makeCntTable(repeats, repeats->htAttributes, "repName");
    cntClass = makeCntTable(repeats, repeats->htAttributes, "repClass");
    cntFamily = makeCntTable(repeats, repeats->htAttributes, "repFamily");

    while(1) {
        pthread_mutex_lock(&fp_lock);
        al = getGroup(al, fp, hdr);
        pthread_mutex_unlock(&fp_lock);
        if(al == NULL) break;
        al_sort(al);
        if(!overlapsGenes(genes, hdr, al, cntTop(al), GTF_IGNORE_STRAND, GTF_MATCH_ANY)) {
            overlapsRepeats(repeats, hdr, al, cntTop(al), GTF_IGNORE_STRAND, GTF_MATCH_ANY, cntName, cntClass, cntFamily);
        }
    }

    //merge cnts
    pthread_mutex_lock(&num_lock);
    mergeCntTable(g_cntName, cntName);
    mergeCntTable(g_cntClass, cntClass);
    mergeCntTable(g_cntFamily, cntFamily);
    pthread_mutex_unlock(&num_lock);

    //Clean up
    destroyCntTable(cntName);
    destroyCntTable(cntClass);
    destroyCntTable(cntFamily);

    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    htsFile *fp;
    bam_hdr_t *hdr;
    GTFtree *genes = NULL, *repeats = NULL;
    aList *al;
    int nThreads = 32;
//    cntTable *cntName, *cntFamily, *cntClass;
    pthread_mutex_init(&fp_lock, NULL);
    pthread_mutex_init(&num_lock, NULL);
    pthread_t threads[nThreads];

    uniqueName = 0, dupName = 0;
    uniqueClass = 0, dupClass = 0;
    uniqueFamily = 0, dupFamily = 0;

    if(argc != 4) {
        fprintf(stderr, "Usage: %s alignments.bam genes.gtf repeats.rmsk\n", argv[0]);
        return 1;
    }

    fp = sam_open(argv[1], "rb");
    hdr = sam_hdr_read(fp);

    //Construct the gene/repeat trees
    genes = GTF2Tree(argv[2], exons);
    fprintf(stderr, "Parsed %s\n", argv[2]); fflush(stderr);
    repeats = RMSK2Tree(argv[3], passThrough);
    fprintf(stderr, "Parsed %s\n", argv[3]); fflush(stderr);
    sortGTF(genes);
    fprintf(stderr, "Converted %s to a tree\n", argv[2]); fflush(stderr);
    sortGTF(repeats);
    fprintf(stderr, "Converted %s to a tree\n", argv[3]); fflush(stderr);

    //Process each chunk of alignments
    al = al_init();
    g_cntName = makeCntTable(repeats, repeats->htAttributes, "repName");
    g_cntClass = makeCntTable(repeats, repeats->htAttributes, "repClass");
    g_cntFamily = makeCntTable(repeats, repeats->htAttributes, "repFamily");
    int i;
    worker_in *win = calloc(1, sizeof(worker_in));
    win->genes = genes;
    win->repeats = repeats;
    win->fp = fp;
    win->hdr = hdr;
    for(i=0; i<nThreads; i++) {
        pthread_create(&threads[i], NULL, worker, (void*) win);
    }
    for(i=0; i<nThreads; i++) {
        pthread_join(threads[i], NULL);
    }

    //Print some tables
    printCntTable(g_cntName, "repName");
    printCntTable(g_cntClass, "repClass");
    printCntTable(g_cntFamily, "repFamily");

    //Free up the cntTables
    destroyCntTable(g_cntName);
    destroyCntTable(g_cntClass);
    destroyCntTable(g_cntFamily);

    //Destroy the trees
    destroyGTFtree(genes);
    destroyGTFtree(repeats);

    //Close things up
    sam_close(fp);
    pthread_mutex_destroy(&fp_lock);
    pthread_mutex_destroy(&num_lock);

    return 0;
}
