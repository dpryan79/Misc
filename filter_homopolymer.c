#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAXREAD 1024

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
    int homopolymer;    //Has the read been homopolymer filtered
} seq_read;

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

void filterHomopolymer(seq_read *read, int maxL) {
    int i = 0, L = 0;
    int l = strlen(read->sequence1);
    char c = NULL;

    read->homopolymer = 0;
    if(l > 0) {
        for(i=0; i<l; i++) {
            if(read->sequence1[i] == c) L++;
            if(L > maxL) {
                read->homopolymer = 1;
                return;
            }
            if(read->sequence1[i] !=c) {
                L = 1;
                c = read->sequence1[i];
            }
        }
    }

    //read #2
    l = strlen(read->sequence2);
    read->homopolymer = 0; 
    if(l > 0) {
        for(i=0; i<l; i++) {
            if(read->sequence2[i] == c) L++;
            if(L > maxL) {
                read->homopolymer = 1;
                return;
            }
            if(read->sequence2[i] !=c) {
                L = 1;
                c = read->sequence2[i];
            }
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
.fastq.gz extension and replacing it with .filtered.fq.gz\n\
\n\
    --homopolymer The maximum allowed length for a homopolymer repeat. If a\n\
             stretch greater than this length is present then the read/pair\n\
             will be filtered. The default is 7.\n\
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
    oname = malloc(sizeof(char)*(strlen(basename) + strlen(".filtered.fq.gz") + 1));
    sprintf(oname, "%s.filtered.fq.gz", basename);

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
    char *line = malloc(MAXREAD*sizeof(char));
    int maxline = MAXREAD;  //This is the buffer size that gzgets will read into
                            //it can be increased internally
    int total_reads = 0, maxL = 7;
    int total_filtered = 0;
    int keep = 0;
    FILE *f1, *f2, *of1, *of2;
    seq_read *read = malloc(sizeof(seq_read));
    char *p = line, *p2;
    char *cmd = malloc(sizeof(char)*MAXREAD);

    //Read in the input
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return(1);
        }
        if(strcmp(argv[i], "--homopolymer") == 0) {
            i++;
            maxL= atoi(argv[i]);
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

        read->homopolymer = 0;
        filterHomopolymer(read, maxL);

        if(read->homopolymer == 1) {
            total_filtered++;
        } else {
            write_trimmed(of1, read, 0);
            if(file2 != NULL) {
                write_trimmed(of2, read, 1);
            }
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
    printf("\t%i\thomopolymer filtered\n", total_filtered);

    //Close things up
    read_destroy(read);
    pclose(f1);
    pclose(of1);
    if(file2 != NULL) {
        pclose(f2);
        pclose(of2);
    }
    free(line);
    free(p);

    return 0;
}
