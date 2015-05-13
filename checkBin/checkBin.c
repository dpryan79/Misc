#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "htslib/sam.h"

//This is directly from the spec.
int calcBin(bam1_t *b) {
    int32_t begin = b->core.pos;
    int32_t end = bam_endpos(b);

    if (begin>>14 == end>>14) return ((1<<15)-1)/7 + (begin>>14);
    if (begin>>17 == end>>17) return ((1<<12)-1)/7 + (begin>>17);
    if (begin>>20 == end>>20) return ((1<<9)-1)/7  + (begin>>20);
    if (begin>>23 == end>>23) return ((1<<6)-1)/7  + (begin>>23);
    if (begin>>26 == end>>26) return ((1<<3)-1)/7  + (begin>>26);
    return 0;
}

void usage(char *prog) {
    fprintf(stderr, "Usage: %s file.bam\n", prog);
    fprintf(stderr, "\n"
"This program accepts a BAM or CRAM file and, for each entry, computes whether\n"
"the associated 'bin' field is correct. If not, the read name and position are\n"
"printed to stdout.\n"
    );
}

int main(int argc, char *argv[]) {
    htsFile *fp;
    bam_hdr_t *hdr;
    bam1_t *b;

    if(argc != 2 || strcmp(argv[1],"-h") == 0) {
        usage(argv[0]);
        return 1;
    }

    fp = hts_open(argv[1], "r");
    if(!fp) {
        fprintf(stderr, "Could not open %s!\n", argv[1]);
        return 1;
    }
    hdr = sam_hdr_read(fp);
    if(!hdr) {
        fprintf(stderr, "Could not read the header from %s!\n", argv[1]);
        return 1;
    }
    b = bam_init1();
    if(!b) {
        fprintf(stderr, "Could not allocate space to hold even a single alignment!\n");
        return 1;
    }
    while(sam_read1(fp, hdr, b) > 0) {
        if(b->core.bin != calcBin(b)) {
            printf("%s @%s:%"PRId32" has bin %i but should be in %i\n", bam_get_qname(b), hdr->target_name[b->core.tid], b->core.pos, b->core.bin, calcBin(b));
        }
    }
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}
