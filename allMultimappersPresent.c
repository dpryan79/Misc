//gcc -O3 -Wall -I/home/ryand/include -L/home/ryand/lib -o /home/ryand/bin/allMultimappersPresent allMultimappersPresent.c -lhts -lz -lpthread
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/hts.h"

typedef struct {
    int l, m;
    bam1_t **b;
} alignmentArray;

//Reset the length
void aaReset(alignmentArray *arr) {
    arr->l = 0;
}

//Increase the size of an alignmentArray
void aaResize(alignmentArray *arr, int l) {
    int i;
    arr->b = realloc(arr->b, sizeof(bam1_t*)*l);
    assert(arr->b);
    for(i=arr->m; i < l; i++) {
        arr->b[i] = bam_init1();
        assert(arr->b[i]);
    }
    arr->m = l;
}

void aaPush(alignmentArray *arr, bam1_t *b) {
    int newM;

    //Reset the size of needed
    if(arr->l == arr->m) {
        newM = arr->l+1;
        kroundup32(newM);
        aaResize(arr, newM);
    }
    bam_copy1(arr->b[arr->l++], b);
}

void aaDestroy(alignmentArray *arr) {
    int i;
    for(i=0; i<arr->m; i++) bam_destroy1(arr->b[i]);
    free(arr->b);
    free(arr);
}

alignmentArray *aaInit() {
    alignmentArray *arr = calloc(1, sizeof(alignmentArray));
    assert(arr);
    return arr;
}

int getNH(bam1_t *b) {
    uint8_t *p = bam_aux_get(b, "NH");
    if(!p) return 1;
    return bam_aux2i(p);
}

int getHI(bam1_t *b) {
    uint8_t *p = bam_aux_get(b, "HI");
    if(!p) return 1;
    return bam_aux2i(p);
}

//Return 1 if at least one alignment is missing, otherwise 0
//If there are no missing alignments, the stack is written to "of"
uint32_t processStack(htsFile *of, bam_hdr_t *hdr, alignmentArray *arr) {
    int NH, i;

    assert(arr->l);
    NH = getNH(arr->b[arr->l-1]);

    if(NH != arr->l && NH>0) return 1;
    for(i=0; i<arr->l; i++) {
        sam_write1(of, hdr, arr->b[i]);
    }
    return 0;
}

void processFile(htsFile *of, bam_hdr_t *hdr, htsFile *fp) {
    alignmentArray *arr = aaInit();
    char *lname = NULL;
    bam1_t *b = bam_init1();
    uint32_t nFiltered = 0, nTotal = 0;

    while(sam_read1(fp, hdr, b) >= 0) {
        if(!lname) {
            aaPush(arr, b);
            lname = bam_get_qname(arr->b[0]);
        } else if(strcmp(lname, bam_get_qname(b)) == 0) {
            aaPush(arr, b);
        } else {
            nTotal++;
            nFiltered += processStack(of, hdr, arr);
            aaReset(arr);
            aaPush(arr, b);
            lname = bam_get_qname(arr->b[0]);
        }
    }
    nTotal++;
    nFiltered += processStack(of, hdr, arr);
    aaDestroy(arr);
    bam_destroy1(b);
fprintf(stderr, "Filtered out %"PRIu32" of %"PRIu32 " groups of multimappers (%5.2f%%)\n", nFiltered, nTotal, 100.0 * ((float) nFiltered)/((float) nTotal));
}

void usage() {
    fprintf(stderr, "Usage: allMultimappersPresent [OPTIONS] <alignments.bam>\n");
    fprintf(stderr, "\n"
"This program will accept a BAM, SAM or CRAM file as input and use the NH and HI\n"
"auxiliary tags to filter out an group of multimapping alignments where at least\n"
"one alignment has been filtered out. This is useful for filtering multimappers\n"
"where at least one alignment is to a feature and has consequently been removed\n"
"(e.g., with bedtools intersect -abam ...). You can use a pipe if the input file\n"
"is '-'. Output is written to stdout, unless otherwise specified.\n"
"\nOPTIONS:\n"
"-o FILE  Output file name (default is to write to stdout)\n"
"-@ INT   Number of compression threads (default: 1)\n"
"-f FILE  File name of the reference genome, for use with CRAM input.\n"
"-B       Output in BAM format, regardless of the input format.\n"
"-C       Output in CRAM format, regardless of the input format (requires -f).\n"
);
}

int main(int argc, char *argv[]) {
    htsFile *fp, *of;
    char *oname = "-", *fa = NULL;
    int nt = 1;
    bam_hdr_t *hdr;
    char c;
    int outType = 0; //1: BAM, 2: CRAM

    opterr = 0;
    while((c = getopt(argc, argv, "hf:@:CBo:")) >= 0) {
        switch(c) {
        case 'h' :
            usage(argv[0]);
            return 0;
        case 'f' :
            fa = optarg;
            break;
        case '@' :
            nt = atoi(optarg);
            if(nt<1) nt=1;
            break;
        case 'C' :
            outType = 2;
            break;
        case 'B' :
            outType = 1;
            break;
        case 'o' :
            oname = optarg;
            break;
        case '?' :
        default :
            fprintf(stderr, "Invalid option '%s'\n", argv[optind-1]);
            usage(argv[0]);
            return 1;
            break;
        }
    }

    if(argc == 1) {
        usage(argv[0]);
        return 0;
    }
    if(argc-optind != 1) {
        fprintf(stderr, "You must specify an input file!\n");
        usage(argv[0]);
        return 1;
    }

    //Open input and output files
    fp = sam_open(argv[optind], "r");
    assert(fp);
    if((outType==2 || fp->is_cram) && !fa) {
        fprintf(stderr, "When either the input or output is in CRAM format, you must specify a reference sequence with -f!\n");
        return 1;
    }
    hdr = sam_hdr_read(fp);
    if(outType) {
        if(outType==2) {
            of = sam_open(oname, "wc");
        } else {
            of = sam_open(oname, "wb");
        }
    } else {
        if(fp->is_cram) {
            of = sam_open(oname, "wc");
        } else {
            of = sam_open(oname, "wb");
        }
    }
    if(fa) {
        hts_set_fai_filename(fp, fa);
        hts_set_fai_filename(of, fa);
    }
    if(nt>1) hts_set_threads(of, nt);
    sam_hdr_write(of, hdr);

    processFile(of, hdr, fp);

    bam_hdr_destroy(hdr);
    sam_close(fp);
    sam_close(of);

    return 0;
}
