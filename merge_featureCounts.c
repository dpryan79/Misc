#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXLINE 1024*1024

int main(int argc, char *argv[]) {
    FILE **fp;
    int i, nfiles = argc-1, finished = 0;
    char *line, *p;

    if(argc<3) {
        fprintf(stderr, "Usage: %s file1.counts file2.counts file3.counts ... > merged.counts\n", argv[0]);
        fprintf(stderr, "\nThis program is aimed at merging count files produced by feature counts into a single matrix.\n");
        return 1;
    }

    fp = malloc(nfiles*sizeof(FILE*));
    for(i=1; i<argc; i++) {
        fp[i-1] = fopen(argv[i], "r");
    }
    line = malloc(sizeof(char)*MAXLINE); //We should really grow this as needed, but this should normally suffice


    //Get past the 2 header lines
    for(i=0; i<nfiles; i++) if(fgets(line, MAXLINE, fp[i]) == NULL) break;
    for(i=0; i<nfiles; i++) if(fgets(line, MAXLINE, fp[i]) == NULL) break;

    //Print a header
    printf("ID");
    for(i=1; i<argc; i++) printf("\t%s", argv[i]);
    printf("\n");

    //Process the important part
    while(1) {
        for(i=0;i<nfiles;i++) {
            if(fgets(line, MAXLINE, fp[i]) == NULL) finished = 1;
            if(finished==1) break;
            p = strtok(line, "\t"); //ID
            //handle the first file differently
            if(i==0) printf("%s", p);
            p = strtok(NULL, "\t"); //Chr
            p = strtok(NULL, "\t"); //Start
            p = strtok(NULL, "\t"); //End
            p = strtok(NULL, "\t"); //Strand
            p = strtok(NULL, "\t"); //Length
            p = strtok(NULL, "\n"); //Count!!!
            if(p==NULL) {
                fprintf(stderr, "Truncated line?\n");
                finished=1;
                break;
            }
            printf("\t%s", p);
        }
        if(finished==1) break;
        printf("\n");
    }

    //Clean things up
    for(i=0;i<nfiles;i++) fclose(fp[i]);
    free(fp);
    free(line);

    return 0;
}
