//gcc -O3 -o bedGraphCorrelogram bedGraphCorrelogram.c -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#define roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

struct call {
    char *chrom;
    int32_t start, stop;
    double meth;
    struct call *next;
};
struct call *first, *last;

typedef struct {
    int l,m;
    double *callsX;
    double *callsY;
} lagArray;

lagArray **lags;

double mean(double *vals, int l) {
    double o = 0.0;
    int i;
    for(i=0; i<l;i++) {
        o += vals[i]/l;
    }
    return o;
}

double stdev(double *vals, int l, double mu) {
    int i;
    double o = 0.0;
    for(i=0; i<l; i++) o += pow(vals[i]-mu, 2.0)/l;
    return sqrt(o);
}

double getNumerator(double *valsX, double *valsY, int l, double muX, double muY) {
    double o = 0.0;
    int i;

    for(i=0; i<l; i++) o += ((valsX[i]-muX)*(valsY[i]-muY))/l;
    return o;
}

void calculateCorsAndOutput(int windowLen) {
    int i;
    double muX, muY, stdevX, stdevY;
    double numerator;
    for(i=0; i<windowLen; i++) {
        if(!lags[i]->l) continue;
        muX = mean(lags[i]->callsX, lags[i]->l);
        muY = mean(lags[i]->callsY, lags[i]->l);
        stdevX = stdev(lags[i]->callsX, lags[i]->l, muX);
        stdevY = stdev(lags[i]->callsY, lags[i]->l, muY);
        numerator = getNumerator(lags[i]->callsX, lags[i]->callsY, lags[i]->l, muX, muY);
        printf("%i\t%f\n", i+1, numerator/(stdevX*stdevY));
    }
    return;
}

lagArray *pushMetrics(lagArray *arr, double a, double b) {
    if(arr->l+1 >= arr->m) {
        arr->m = arr->l+1;
        arr->m = roundup32(arr->m);
        arr->callsX = realloc(arr->callsX, sizeof(double)*arr->m);
        assert(arr->callsX);
        arr->callsY = realloc(arr->callsY, sizeof(double)*arr->m);
        assert(arr->callsY);
    }
    arr->callsX[arr->l] = a;
    arr->callsY[arr->l++] = b;
    return arr;
}

inline int getDistance(struct call *a, struct call *b) {
    return b->start - a->stop + 1;
}

void processLL(int windowLen) {
    struct call *tmp = first->next;
    int dist;
    while(tmp) {
        if(strcmp(first->chrom, tmp->chrom) == 0) {
            dist = getDistance(first, tmp);
            if(dist <= windowLen) {
                lags[dist-1] = pushMetrics(lags[dist-1], first->meth, tmp->meth);
                tmp = tmp->next;
            } else {
                break;
            }
        } else {
            break;
        }
    }
}

int processLine(char *line, struct call *mcall) {
    int nMeth, nUmeth;
    char *p;

    //Chrom
    p = strtok(line, "\t");
    mcall->chrom = strdup(p);
    //Start
    p = strtok(NULL, "\t");
    mcall->start = atoi(p);
    //end
    p = strtok(NULL, "\t");
    mcall->stop = atoi(p);
    //score
    p = strtok(NULL, "\t");
    //nMeth
    p = strtok(NULL, "\t");
    nMeth = atoi(p);
    //nUmeth
    p = strtok(NULL, "\n");
    nUmeth = atoi(p);

    if(nMeth+nUmeth) {
        mcall->meth = ((double) nMeth)/((double)(nMeth+nUmeth));
    } else {
        mcall->meth = 0.0;
    }
    return nMeth+nUmeth;
}

void pushCall(struct call *mcall) {
    if(first == NULL) {
        first = mcall;
        last = mcall;
    } else {
        last->next = mcall;
        last = mcall;
    }
}

void process(FILE *fp, int minDepth, int windowLen) {
    char *line = malloc(sizeof(char)*1024);
    int llLen = 0; //length of the linked-list
    struct call *mcall, *tmp;
    assert(line);

    while(fgets(line, 1024, fp)) {
        if(strncmp("track", line, 5) == 0) continue;

        mcall = calloc(1, sizeof(struct call));
        assert(mcall);
        if(processLine(line, mcall) < minDepth) {
            free(mcall->chrom);
            free(mcall);
            continue;
        }

        pushCall(mcall);
        if(llLen >= windowLen) {
            processLL(windowLen);
            tmp = first;
            first = first->next;
            free(tmp->chrom);
            free(tmp);
        } else {
            llLen++;
        }
    }
    while(first != last) {
        processLL(windowLen);
        tmp = first;
        first = first->next;
        free(tmp->chrom);
        free(tmp);
    }
    free(first->chrom);
    free(first);
    free(line);
}

void usage(char *prog) {
    fprintf(stderr, "Usage: %s [OPTIONS] <file.bedGraph>\n", prog);
    fprintf(stderr, "\n\
%s will take a sorted BedGraph file containing with the number of\n\
methylated Cs in column 5 and the number of unmethylated Cs in column 6. This is\n\
the format output by bison, bismark and PileOMeth.\n", prog);
    fprintf(stderr, "\n\
-d INT	The minimum depth for inclusion (default, 5). A value of 0 will include\n\
	all sites.\n\
-l INT	The maximum distance to compute the correlation over (default, 1000 bp).\n");
}

int main(int argc, char *argv[]) {
    FILE *fp = NULL;
    int c, i, minDepth = 5, windowLen = 1000;

    first = last = NULL;

    while((c = getopt(argc, argv, "hd:l:")) >= 0) {
        switch(c) {
        case 'h' :
            usage(argv[0]);
            return 0;
        case 'd' :
            minDepth = atoi(optarg);
            break;
        case 'l' :
            windowLen = atoi(optarg);
            break;
        default :
            fprintf(stderr, "Invalid option %c\n", c);
            usage(argv[0]);
            return 1;
        }
    }
    if(argc==1) {
        usage(argv[0]);
        return 0;
    }
    if(argc-optind != 1) {
        fprintf(stderr, "Missing input file!\n");
        usage(argv[0]);
        return 1;
    }

    fp = fopen(argv[optind], "r");
    assert(fp);

    lags = calloc(windowLen, sizeof(lagArray *));
    assert(lags);
    for(i=0; i<windowLen; i++) {
        lags[i] = calloc(1, sizeof(lagArray));
        assert(lags[i]);
    }

    process(fp, minDepth, windowLen);
    calculateCorsAndOutput(windowLen);

    for(i=0; i<windowLen; i++) {
        if(lags[i]->l) {
            free(lags[i]->callsX);
            free(lags[i]->callsY);
        }
        free(lags[i]);
    }
    free(lags);
    fclose(fp);

    return 0;
}
