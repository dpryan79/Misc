#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAXREAD 1024
#define DELETION_COST 1
#define INSERTION_COST 1
#define MISMATCH_COST 1
#define MATCH_COST 0
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

//Trim from the right: 6,7, 10, 11, 14,15, 0->start2
//Trim from the left: 10-15, stop2->strlen(read)
enum FLAG {
    START_WITHIN_SEQ1 = 1,
    START_WITHIN_SEQ2 = 2,
    STOP_WITHIN_SEQ1 = 4,
    STOP_WITHIN_SEQ2 = 8,
    SEMIGLOBAL = 15,
    ALLOW_WILDCARD_SEQ1 = 1,
    ALLOW_WILDCARD_SEQ2 = 2,
};

typedef struct {
    int flags;
    float error_rate;
    int min_overlap;
    int min_length;
    int keep;
    int base;
    int cutoff;
    int rrbs;
    int mspi;
    int taqi;
    int non_directional;
    int polyA;
} config;

typedef struct {
    char *name1;	//Name of left read, max length of maxname1
    char *sequence1;	//Sequence of the left read, max length of maxread1
    char *qual1;	//Quality score of the left read, max length of maxread1
    char *name2;	//As above, but for the left read
    char *sequence2;
    char *qual2;
    int maxname1;
    int maxname2;
    int maxread1;
    int maxread2;
    int quality_trimmed;	//Has the read been quality trimmed? 0: No, 1: Yes
    int adapter_trimmed;	//Has the read been adapter trimmed? 0: No, 1: Yes
    int Threeprime_adapter_trimmed1;	//This is only needed for RRBS!
    int Threeprime_adapter_trimmed2;	//This is only needed for RRBS!
    int Fiveprime_adapter_trimmed1;	//This could be needed in the future
    int Fiveprime_adapter_trimmed2;	//This could be needed in the future
    int polyA_trimmed;		//Has the read been 3' polyA trimmed? 0: No, 1: Yes
} seq_read;

//See cutadapt
typedef struct {
    int cost;
    int matches;
    int origin;
} Entry;

typedef struct {
    int start1;
    int stop1;
    int start2;
    int stop2;
    int matches;
    int errors;
} alignment;

//Allocate space required by a read
//return 0 on error
int read_init(seq_read *read) {
    int rv = 1;
    read->name1 = calloc(sizeof(char), MAXREAD);
    read->name2 = calloc(sizeof(char), MAXREAD);
    read->sequence1 = calloc(sizeof(char), MAXREAD);
    read->sequence2 = calloc(sizeof(char), MAXREAD);
    read->qual1 = calloc(sizeof(char), MAXREAD);
    read->qual2 = calloc(sizeof(char), MAXREAD);

    read->maxname1 = MAXREAD;
    read->maxname2 = MAXREAD;
    read->maxread1 = MAXREAD;
    read->maxread2 = MAXREAD;

    if(read->name1 == NULL) rv=0;
    if(read->name2 == NULL) rv=0;
    if(read->sequence1 == NULL) rv=0;
    if(read->sequence2 == NULL) rv=0;
    if(read->qual1 == NULL) rv=0;
    if(read->qual2 == NULL) rv=0;
    return rv;
}

//Free up the memory occupied by a read
void read_destroy(seq_read *read) {
    free(read->name1);
    free(read->name2);
    free(read->sequence1);
    free(read->sequence2);
    free(read->qual1);
    free(read->qual2);
    free(read);
}

//Reverse complement an adapter or other sequence
//The output must be free()d
char *reverse_complement(char *seq) {
    int i;
    int maxi = strlen(seq);
    char c;
    char *output = calloc(sizeof(char), maxi+1);

    for(i=maxi; i > 0; i--) {
        switch(seq[i-1]) {
          case 'A' :
            c = 'T';
            break;
          case 'C' :
            c = 'G';
            break;
          case 'G' :
            c = 'C';
            break;
          case 'T' :
            c = 'A';
            break;
          case 'N' :
            c = 'N';
            break;
          default :
            c = 'N';
            printf("An unexpected nucleotide, %c, was encountered in a read. It was replaced with an N.\n", c);
        }
        output[maxi-i] = c;
    }
    return output;
}

//This is effectively copied from calignmodule.c in cutadapt
//the output needs to be free()d
alignment* global_alignment(char *s1, char *s2, config config) {
    int i, j, best_i, best_j, best_cost, best_matches, best_origin;
    int k, last;
    int m = strlen(s1), n = strlen(s2);
    int degenerate = 0;
    int match, cost_diag, cost_deletion, cost_insertion, origin, matches;
    int length, cost;
    int start1, start2;
    alignment *output;
    Entry *column, tmp_entry;

    /*
    DP Matrix:
              s2 (j)
            ----------> n
           |
    s1 (i) |
           |
           V
           m
    */

    // only a single column of the DP matrix is stored
    column = (Entry*)malloc((m+1)*sizeof(Entry));
    if (column == NULL) return NULL;
    output = malloc(sizeof(alignment));

    //initialize first column
    for (i = 0; i <= m; ++i) {
        column[i].matches = 0;
        column[i].cost = (config.flags & START_WITHIN_SEQ1) ? 0 : i * DELETION_COST;
        column[i].origin = (config.flags & START_WITHIN_SEQ1) ? -i : 0;
    }

    best_i = m;
    best_j = 0;
    best_cost = column[m].cost;
    best_matches = 0;
    best_origin = column[m].origin;

    // maximum no. of errors
    k = config.error_rate * m;
    last = k + 1;
    if (config.flags & START_WITHIN_SEQ1) {
        last = m;
    }
    // iterate over columns
    for (j = 1; j <= n; ++j) {
        // remember first entry
        tmp_entry = column[0];

        // fill in first entry in this column TODO move out of loop
        if (config.flags & START_WITHIN_SEQ2) {
            column[0].cost = 0;
            column[0].origin = j;
            column[0].matches = 0;
        } else {
            column[0].cost = j * INSERTION_COST;
            column[0].origin = 0;
            column[0].matches = 0;
        }
        for (i = 1; i <= last; ++i) {
            match = (s1[i-1] == s2[j-1])
                || ((degenerate & ALLOW_WILDCARD_SEQ1) && (s1[i-1] == 'N'))
                || ((degenerate & ALLOW_WILDCARD_SEQ2) && (s2[j-1] == 'N'));
            cost_diag = tmp_entry.cost + (match ? MATCH_COST : MISMATCH_COST);
            cost_deletion = column[i].cost + DELETION_COST;
            cost_insertion = column[i-1].cost + INSERTION_COST;

            if (cost_diag <= cost_deletion && cost_diag <= cost_insertion) {
                // MATCH or MISMATCH
                cost = cost_diag;
                origin = tmp_entry.origin;
                matches = tmp_entry.matches + match;
            } else if (cost_insertion <= cost_deletion) {
                // INSERTION
                cost = cost_insertion;
                origin = column[i-1].origin;
                matches = column[i-1].matches;
            } else {
                // DELETION
                cost = cost_deletion;
                origin = column[i].origin;
                matches = column[i].matches;
            }

            // remember current cell for next iteration
            tmp_entry = column[i];

            column[i].cost = cost;
            column[i].origin = origin;
            column[i].matches = matches;

        }
        while (column[last].cost > k) {
            last--;
        }
        if (last < m) {
            last++;
        } else {
            // found
            // if requested, find best match in last row
            if (config.flags & STOP_WITHIN_SEQ2) {
                // length of the aligned part of string1
                length = m + min(column[m].origin, 0);
                cost = column[m].cost;
                matches = column[m].matches;
                if (cost <= length * config.error_rate && (matches > best_matches || (matches == best_matches && cost < best_cost))) {
                    // update
                    best_matches = matches;
                    best_cost = cost;
                    best_origin = column[m].origin;
                    best_i = m;
                    best_j = j;
                }
            }

        }
        // column finished
    }

    if (config.flags & STOP_WITHIN_SEQ1) {
        // search in last column // TODO last?
        for (i = 0; i <= m; ++i) {
            length = i + min(column[i].origin, 0);
            cost = column[i].cost;
            matches = column[i].matches;
            if (cost <= length * config.error_rate && (matches > best_matches || (matches == best_matches && cost < best_cost))) {
                // update best
                best_matches = matches;
                best_cost = cost;
                best_origin = column[i].origin;
                best_i = i;
                best_j = n;
            }
        }
    }

    free(column);
    if (best_origin >= 0) {
        start1 = 0;
        start2 = best_origin;
    } else {
        start1 = -best_origin;
        start2 = 0;
    }

    output->start1 = start1;
    output->stop1 = best_i;
    output->start2 = start2;
    output->stop2 = best_j;
    output->matches = best_matches;
    output->errors = best_cost;
    return output;
}

//Need to handle, MspI, TaqI and non_directional
void rrbs_trim(seq_read *read, config config, int side) {
    int len, i, j;

    if(side == 3) {
        //Currently, we only have methods for MspI and/or TaqI
        //Digested libraries
        if(config.mspi || config.taqi) {
            if(read->Threeprime_adapter_trimmed1) {
                len = strlen(read->sequence1);
                if(len >= 2) {
                    read->sequence1[len-2] = '\0';
                    read->qual1[len-2] = '\0';
                }
            }
            if(read->Threeprime_adapter_trimmed2) {
                len = strlen(read->sequence2);
                if(len >= 2) {
                    read->sequence2[len-2] = '\0';
                    read->qual2[len-2] = '\0';
                }
            }
        }
    } else if(side == 5) {
        //After the first round of quality and adapter trimming, trim off 2bp if the read starts with CGA or CAA
        if(config.non_directional && (config.mspi || config.taqi)) {
            if(strncmp(read->sequence1, "CGA", 3) == 0 || strncmp(read->sequence1, "CAA", 3) == 0) {
                len = strlen(read->sequence1);
                for(i=0, j=2; j<len; i++, j++) {
                    read->sequence1[i] = read->sequence1[j];
                    read->qual1[i] = read->qual1[j];
                }
                read->sequence1[i] = '\0';
                read->qual1[i] = '\0';
            }
            if(strncmp(read->sequence2, "CGA", 3) == 0 || strncmp(read->sequence2, "CAA", 3) == 0) {
                len = strlen(read->sequence2);
                for(i=0, j=2; j<len; i++, j++) {
                    read->sequence2[i] = read->sequence2[j];
                    read->qual2[i] = read->qual2[j];
                }
                read->sequence2[i] = '\0';
                read->qual2[i] = '\0';
            }
        }
    }
}
        

//Match an adapter and return the sequence, assuming the adapter is on the 3' end
//Return 1 if trimmed and 0 if not
int trim_3prime(char *adapter, seq_read *read, config config) {
    alignment *align1, *align2;
    char *adapter2 = reverse_complement(adapter);
    int rv = 0;

    align1 = global_alignment(adapter2, read->sequence1, config);
    if(align1->matches >= config.min_overlap) {
        read->sequence1[align1->start2] = '\0';
        read->qual1[align1->start2] = '\0';
        read->Threeprime_adapter_trimmed1 = 1;
        rv = 1;
    }
    free(align1);

    if(strlen(read->sequence2)) {
        align2 = global_alignment(adapter2, read->sequence2, config);
        if(align2->matches >= config.min_overlap) {
            read->sequence2[align2->start2] = '\0';
            read->qual2[align2->start2] = '\0';
            read->Threeprime_adapter_trimmed2 = 1;
            rv = 1;
        }
        free(align2);
    }

    free(adapter2);
    return rv;
}

int trim_5prime(char *adapter, seq_read *read, config config) {
    alignment *align1, *align2;
    int i, j, rv = 0;

    align1 = global_alignment(adapter, read->sequence1, config);
    if(align1->matches >= config.min_overlap) {
        for(i=align1->stop2, j=0; i<strlen(read->sequence1); i++, j++) {
            read->sequence1[j] = read->sequence1[i];
            read->qual1[j] = read->qual1[i];
        }
        read->sequence1[j] = '\0';
        read->qual1[j] = '\0';
        read->Fiveprime_adapter_trimmed1 = 1;
        rv = 1;
    }
    free(align1);

    //If this is a paired-end read, do the same
    if(strlen(read->sequence2)) {
        align2 = global_alignment(adapter, read->sequence2, config);
        if(align2->matches >= config.min_overlap) {
            for(i=align2->stop2, j=0; i<strlen(read->sequence2); i++, j++) {
                read->sequence2[j] = read->sequence2[i];
                read->qual2[j] = read->qual2[i];
            }
            read->sequence2[j] = '\0';
            read->qual2[j] = '\0';
            read->Fiveprime_adapter_trimmed2 = 1;
            rv = 1;
        }
        free(align2);
    }

    return rv;
}

void quality_trim(seq_read *read, config config) {
    int i, j, cutoff;
    int trimmed = 0;

    //First, go through the 3' end
    cutoff=strlen(read->qual1);
    for(i=strlen(read->qual1)-1; i>0; i--) {
        if(read->sequence1[i] == 'N') {
            cutoff--;
        } else if((read->qual1[i] - config.base) <= config.cutoff) {
            cutoff--;
        } else if((read->qual1[i] + read->qual1[i-1] - 2*config.base) <= 2*config.cutoff) {
            cutoff--;
        } else {
            break;
        }
    }
    if(cutoff != strlen(read->sequence1)) trimmed = 1;
    read->sequence1[cutoff] = '\0';
    read->qual1[cutoff] = '\0';
    //5'
    cutoff = 0;
    for(i=0; i<strlen(read->qual1)-1; i++) {
        if(read->sequence1[i] == 'N') {
            cutoff++;
        } else if((read->qual1[i] - config.base) <= config.cutoff) {
            cutoff++;
        } else if((read->qual1[i] + read->qual1[i+1] - 2*config.base) <= 2*config.cutoff) {
            cutoff++;
        } else {
            break;
        }
    }
    if(cutoff > 0) {
        for(i=cutoff, j=0; i<strlen(read->qual1); i++, j++) {
            read->qual1[j] = read->qual1[i];
            read->sequence1[j] = read->sequence1[i];
        }
        read->qual1[j] = '\0';
        read->sequence1[j] = '\0';
        trimmed = 1;
    }

    //Is this paired-ended?
    if(strlen(read->name2)) {
        cutoff=strlen(read->qual2);
        for(i=strlen(read->qual2)-1; i>0; i--) {
            if(read->sequence2[i] == 'N') {
                cutoff--;
            } else if((read->qual2[i] - config.base) <= config.cutoff) {
                cutoff--;
            } else if((read->qual2[i] + read->qual2[i-1] - 2*config.base) <= 2*config.cutoff) {
                cutoff--;
            } else {
                break;
            }
        }
        if(cutoff != strlen(read->sequence1)) trimmed = 1;
        read->sequence2[cutoff] = '\0';
        read->qual2[cutoff] = '\0';
        if(cutoff > 0) { //For some reason, this segfaults if the read is empty
            cutoff = 0;
            for(i=0; i<strlen(read->qual2)-1; i++) {
                if(read->sequence2[i] == 'N') {
                    cutoff++;
                } else if((read->qual2[i] - config.base) <= config.cutoff) {
                    cutoff++;
                } else if((read->qual2[i] + read->qual2[i+1] - 2*config.base) <= 2*config.cutoff) {
                    cutoff++;
                } else {
                    break;
                }
            }
        }
        if(cutoff > 0) {
            trimmed = 1;
            for(i=cutoff, j=0; i<strlen(read->qual2); i++, j++) {
                read->qual2[j] = read->qual2[i];
                read->sequence2[j] = read->sequence2[i];
            }
            read->qual2[j] = '\0';
            read->sequence2[j] = '\0';
        }
    }

    if(trimmed) read->quality_trimmed = trimmed;
}

void trimPolyA(seq_read *read, config *config) {
    int l = strlen(read->sequence1), i, j = 0;

    read->polyA_trimmed = 0;
    if(l > 0) {
        j = 0;
        for(i = l-1; i > 0; i--) {
            if(read->sequence1[i] == 'A') j++;
            else break;
        }
        if(j > config->polyA) {
            read->sequence1[l-j] = '\0';
            read->qual1[l-j] = '\0';
            read->polyA_trimmed = 1;
        }
    }

    l = strlen(read->sequence2);
    if(l > 0) {
        j = 0;
        for(i = l-1; i > 0; i--) {
            if(read->sequence2[i] == 'A') j++;
            else break;
        }
        if(j > config->polyA) {
            read->sequence2[l-j] = '\0';
            read->qual2[l-j] = '\0';
            read->polyA_trimmed = 1;
        }
    }
}

void usage(char *prog) {
    printf("%s [options] [-1 fastq_1.gz -2 fastq_2.gz | fastq.gz]\n", prog);
    printf("\n\
Input files are gzipped fastq files. If reads are paired-end rather than\n\
single-ended, -1 and -2 must be used to denote the left and right reads\n\
(usually ending in _1.fq.gz and _2.fq.gz or similar, respectively)\n\
\n\
Output files will be gzipped and named by removing the .fq.gz or \n\
.fastq.gz extension and replacing it with .trimmed.fq.gz\n\
\n\
    -e       Adapter alignment error rate, default is 0.15\n\
\n\
    -a       Adapter sequence. N.B., the sequence is as appended to the\n\
             5' end of the left read of a fragment! Many do the opposite\n\
             (i.e., they use the reverse complement). The default is the\n\
             Illumina adapter, GCTCTTCCGATCT\n\
\n\
    -overlap Minimum overlap of the adapter with a read to result in\n\
             trimming. Default is 3. Lower values are likely only needed\n\
             in RRBS\n\
\n\
    -min_length  Minumum read length for output. If either read of a\n\
             paired-end set is < min_length, then neither will be output\n\
             (see -keep). Default is 20.\n\
\n\
    --keep   Keep reads whose mate was too short for output (writtend to\n\
             .orphaned.fq.gz files\n\
\n\
    -q       Quality score threshold. Default is 20.\n\
\n\
    --64     Quality scores are Phred+64 rather than Phred+33 (default)\n\
\n\
    --MspI   This is an RRBS library digested with MspI (possible in\n\
             addition to other enzymes). Reads in which an adapter is\n\
             trimmed from the 3' side will have an additional 2bp\n\
             trimmed from their 3' sides.\n\
\n\
    --TaqI   As with --MspI, these can both be specified and are both\n\
             currently treated the same.\n\
\n\
    --non_directional  This is a non-directional RRBS library.\n\
             Currently, --MspI and/or --TaqI must be specified for this\n\
             to have an effect. After quality and 5' adapter trimming,\n\
             reads beginning with CGA or CAA (i.e., possibly starting\n\
             with unreliable filled in bases) will have the first 2 bp\n\
             removed.\n\
\n\
    --polyA  The maximum number of 3' As at the end of a read to allow. If more\n\
             than this are present, then all As will be trimmed off. A value of\n\
             0 (the default) indicates to ignore polyAs.\n\
\n");
}

FILE* determine_name(char *fname, int orphaned) {
    char *basename = malloc(sizeof(char)*MAXREAD);
    char *p = basename;
    char *oname;
    char *cmd = malloc(sizeof(char)*MAXREAD);
    FILE *output;

    //This is unlikely to ever occur
    if(strlen(fname) > MAXREAD-1) {
        basename = realloc(basename, sizeof(char) * (strlen(fname) + 20));
        printf("determine_name\n");
        if(p!=basename) free(p);
    }
    basename = strcpy(basename, fname);

    p = strrchr(basename, '.');
    if(p != NULL) {
        if(strcmp(p, ".gz") == 0 || strcmp(p, ".bz2") == 0 || strcmp(p, ".bz") == 0 || strcmp(p, ".GZ") == 0) {
            *p = '\0';
            p = strrchr(basename, '.');
            if(p != NULL) {
                if(strcmp(p, ".fastq") == 0 || strcmp(p, ".fq") == 0) *p = '\0';
            }
        }
    }

    //construct the output file name
    if(orphaned) {
        oname = malloc(sizeof(char)*(strlen(basename) + strlen(".orphaned.fq.gz") + 1));
        sprintf(oname, "%s.orphaned.fq.gz", basename);
    } else {
        oname = malloc(sizeof(char)*(strlen(basename) + strlen(".trimmed.fq.gz") + 1));
        sprintf(oname, "%s.trimmed.fq.gz", basename);
    }

    //Finally, open the file for writing
    sprintf(cmd, "gzip > %s", oname);
    output = popen(cmd, "w");
    free(cmd);

    free(basename);
    free(oname);
    return output;
}

void write_trimmed(FILE *f, seq_read *read, int lr) {
    if(lr == 0) {
        fprintf(f, "%s\n", read->name1);
        fprintf(f, "%s\n", read->sequence1);
        fprintf(f, "+\n");
        fprintf(f, "%s\n", read->qual1);
    } else {
        fprintf(f, "%s\n", read->name2);
        fprintf(f, "%s\n", read->sequence2);
        fprintf(f, "+\n");
        fprintf(f, "%s\n", read->qual2);
    }
}

int main(int argc, char *argv[]) {
    int i;
    char *file1 = NULL, *file2 = NULL;
    char *adapter = NULL;
    char *line = malloc(MAXREAD*sizeof(char));
    int maxline = MAXREAD;  //This is the buffer size that gzgets will read into
                            //it can be increased internally
    int total_reads = 0;
    int total_adapter_trimmed = 0, total_discarded = 0;
    int total_quality_trimmed = 0, total_orphaned = 0;
    int total_polyA_trimmed = 0;
    int keep = 0;
    FILE *f1, *f2, *of1, *of2, *of1_orphaned = NULL, *of2_orphaned = NULL;
    config config;
    seq_read *read = malloc(sizeof(seq_read));
    char *p = line, *p2;
    char *cmd = malloc(sizeof(char)*MAXREAD);

    //Defaults
    config.error_rate = 0.15;
    config.min_overlap = 3;
    config.flags = 10; //This is not currently user alterable
    config.min_length = 20;
    config.keep = 0;
    config.base = 33;
    config.cutoff = 20;
    config.rrbs = 0;
    config.mspi = 0;
    config.taqi = 0;
    config.non_directional = 0;
    config.polyA = 0;

    //Read in the input
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return(1);
        }
        if(strcmp(argv[i], "-e") == 0) {
            i++;
            config.error_rate = atof(argv[i]);
        } else if(strcmp(argv[i], "--keep") == 0) {
            config.keep = 1;
        } else if(strcmp(argv[i], "--64") == 0) {
            config.base = 64;
        } else if(strcmp(argv[i], "-q") == 0) {
            i++;
            config.cutoff = atoi(argv[i]);
        } else if(strcmp(argv[i], "--MspI") == 0) {
            config.rrbs = 1;
            config.mspi = 1;
        } else if(strcmp(argv[i], "--TaqI") == 0) {
            config.rrbs = 1;
            config.taqi = 1;
        } else if(strcmp(argv[i], "--non_directional") == 0) {
            config.rrbs = 1;
            config.non_directional = 1;
        } else if(strcmp(argv[i], "-min_length") == 0) {
            i++;
            config.min_length = atoi(argv[i]);
        } else if(strcmp(argv[i], "-a") == 0) {
            i++;
            adapter = malloc(sizeof(char) * (strlen(argv[i])+1));
            adapter = strcpy(adapter, argv[i]);
        } else if(strcmp(argv[i], "-overlap") == 0) {
            i++;
            config.min_overlap = atoi(argv[i]);
        } else if(strcmp(argv[i], "-2") == 0) {
            i++;
            file2 = argv[i];
        } else if(strcmp(argv[i], "-1") == 0) {
            i++;
            file1 = argv[i];
        } else if(strcmp(argv[i], "--polyA") == 0) {
            i++;
            config.polyA = atoi(argv[i]);
        } else {
            if(strncmp(argv[i], "-", 1) == 0) {
                //Got an unmatched option
                usage(argv[0]);
                return(1);
            }
            if(file1 == NULL) {
                file1 = argv[i];
            } else {
                printf("To make things simpler, you must specify paired-end reads with -1 and -2.\n");
                usage(argv[0]);
                return(1);
            }
        }
    }
    if(file1 == NULL) {
        usage(argv[0]);
        return 1;
    }

    if(adapter == NULL) {
        adapter = malloc(sizeof(char) * (strlen("GCTCTTCCGATCT")+1));
        adapter = strcpy(adapter, "GCTCTTCCGATCT");
    }

    //Open the files
    p2 = strrchr(file1, '.');
    if(strcmp(p2, ".gz") == 0 || strcmp(p2, ".GZ") == 0) {
        sprintf(cmd, "zcat %s", file1);
    } else if(strcmp(p2, ".bz2") == 0 || strcmp(p2, ".bz") == 0) {
        sprintf(cmd, "bzcat %s", file1);
    } else {
        sprintf(cmd, "cat %s", file1);
    }
    f1 = popen(cmd, "r");
    of1 = determine_name(file1, 0);
    if(keep) of1_orphaned = determine_name(file1, 1);
    if(file2 != NULL) {
        if(strcmp(p2, ".gz") == 0 || strcmp(p2, ".GZ") == 0) {
            sprintf(cmd, "zcat %s", file2);
        } else if(strcmp(p2, ".bz2") == 0 || strcmp(p2, ".bz") == 0) {
            sprintf(cmd, "bzcat %s", file2);
        } else {
            sprintf(cmd, "cat %s", file2);
        }
        f2 = popen(cmd, "r");
        of2 = determine_name(file2, 0);
        if(keep) of2_orphaned = determine_name(file2, 1);
    }
    free(cmd);

    /******************************************************************
    /
    /   Everything Below here can be put in a thread function with
    /   multiple mutexes
    /
    ******************************************************************/
    //Initialize the read
    read_init(read);

    //Parse each read
    while(1) {
        //Lock a mutex
        //Read1 Name
        line = fgets(line, maxline, f1);
        if(line == NULL) break;
        line[strlen(line)-1] = '\0';
        read->name1 = strcpy(read->name1, line);
        //Read1 sequence
        line = fgets(line, maxline, f1);
        line[strlen(line)-1] = '\0';
        read->sequence1 = strcpy(read->sequence1, line);
        //+
        line = fgets(line, maxline, f1);
        line[strlen(line)-1] = '\0';
        //Read1 Quality
        line = fgets(line, maxline, f1);
        line[strlen(line)-1] = '\0';
        read->qual1 = strcpy(read->qual1, line);

        //Read2
        if(file2 != NULL) {
            //Read1 Name
            line = fgets(line, maxline, f2);
            line[strlen(line)-1] = '\0';
            read->name2 = strcpy(read->name2, line);
            //Read1 sequence
            line = fgets(line, maxline, f2);
            line[strlen(line)-1] = '\0';
            read->sequence2 = strcpy(read->sequence2, line);
            //+
            line = fgets(line, maxline, f2);
            line[strlen(line)-1] = '\0';
            //Read1 Quality
            line = fgets(line, maxline, f2);
            line[strlen(line)-1] = '\0';
            read->qual2 = strcpy(read->qual2, line);
        }
        //Unlock a mutex

        read->quality_trimmed = 0;
        read->adapter_trimmed = 0;
        quality_trim(read, config);
        if(trim_5prime(adapter, read, config)) {
            //RRBS? We need to run before quality_trim regardless of whether
            //and adapter was 5' trimmed
            if(config.rrbs) rrbs_trim(read, config, 5);
            //Perform additional quality and N trimming!
            quality_trim(read, config);
            read->adapter_trimmed = 1;
        } else if(config.rrbs) {
            rrbs_trim(read, config, 5);
        }
        if(trim_3prime(adapter, read, config)) {
            //Perform additional quality and N trimming!
            if(config.rrbs) rrbs_trim(read, config, 3);
            quality_trim(read, config);
            read->adapter_trimmed = 1;
        }

        if(config.polyA > 0) {
            trimPolyA(read, &config);
        }

        //Multithreading is done, only the main process should write to files!
        if(read->quality_trimmed) total_quality_trimmed++;
        if(read->adapter_trimmed) total_adapter_trimmed++;
        if(read->polyA_trimmed) total_polyA_trimmed++;
        if(strlen(read->sequence1) >= config.min_length) {
            if(file2 != NULL) { //Paired-end
                if(strlen(read->sequence2) >= config.min_length) {
                    //Lock a mutex
                    write_trimmed(of1, read, 0);
                    //Unlock a mutex
                    //Lock a mutex
                    write_trimmed(of2, read, 1);
                    //Unlock a mutex
                } else if(config.keep) {
                    //Lock a mutex
                    write_trimmed(of1_orphaned, read, 0);
                    //Unlock a mutex
                    total_orphaned++;
                }
            } else {
                //Lock a mutex
                write_trimmed(of1, read, 0);
                //Unlock a mutex
            }
        } else if(config.keep && file2 != NULL) {
            if(strlen(read->sequence2) >= config.min_length) {
                //Lock a mutex
                write_trimmed(of2_orphaned, read, 1);
                //Unlock a mutex
                total_orphaned++;
            } else {
                total_discarded++;
            }
        } else {
            total_discarded++;
        }
        total_reads++;

        //Give some output
        if(total_reads % 1000000 == 0) printf("%i reads processed\n", total_reads);
    }

    /******************************************************************
    /
    /   This marks the end of the threading function
    /
    ******************************************************************/

    //Write some diagnostic output
    printf("There were %i reads processed.\n", total_reads);
    printf("\t%i\tadapter trimmed\n", total_adapter_trimmed);
    printf("\t%i\tquality trimmed\n", total_quality_trimmed);
    printf("\t%i\t3' polyA trimmed\n", total_polyA_trimmed);
    printf("\t%i\torphaned\n", total_orphaned);
    printf("\t%i\tdiscarded\n", total_discarded);

    //Close things up
    read_destroy(read);
    pclose(f1);
    pclose(of1);
    if(file2 != NULL) {
        pclose(f2);
        pclose(of2);
    }
    if(config.keep) {
        pclose(of1_orphaned);
        pclose(of2_orphaned);
    }
    free(line);
    free(p);
    free(adapter);
    
    return 0;
}
