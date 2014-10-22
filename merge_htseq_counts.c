#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#define MAXLINE 512

void print_labels(char *labels) {
    char *label = NULL;

    printf("ID");
    label = strtok(labels,",");
    printf("\t%s",label);
    while(1) {
        label = strtok(NULL,",");
        if(label == NULL) break;
        printf("\t%s",label);
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    FILE **fp_array = NULL;
    char **line_array = NULL;
    char *label = NULL;
    int num_files = 0, i = 0, j = 0, end = 0, count = 0;
    struct stat info;

    /* read in the file names */
    if(argc < 3) {
        fprintf(stderr,"Usage: %s [-L Column,Labels] file1 file2 file3 ...\n \
		\t-L\t(Optional) Comma separated list of column labels. Should be one per file.\n \
		Output is to stdout, so you'll want to redirect the output to a file.\n",argv[0]);
        return 1;
    };

    for(i=1,j=0;i<argc;i++) {
        if(strcmp(argv[i],"-L") == 0) {
            print_labels(argv[i+1]);
            i++;
            continue;
        }
        num_files++;
        fp_array = realloc(fp_array, num_files * sizeof(FILE *));
        line_array = realloc(line_array, num_files * sizeof(char *));
        //Ensure that the file actually exists
        if(stat(argv[i], &info) == 0) {
            fp_array[num_files-1] = fopen(argv[i],"r");
        } else {
            fprintf(stderr,"%s does not exist.\n",argv[i]);
            return 2;
        }
        line_array[num_files-1] = malloc(sizeof(char[MAXLINE]));
    }

    while(1) {
        for(i=0;i<num_files;i++) {
            if(fgets(line_array[i],MAXLINE,fp_array[i]) == NULL) {
                end = 1;
            }
        }
        if(end == 1) break;

        //The first column in each file must match, otherwise this won't work as efficiently
        label = strtok(line_array[0],"\t");
        if(strcmp("no_feature",label) == 0) break;
        count = atoi(strtok(NULL,"\t"));
        printf("%s\t%d",label,count);
        for(i=1;i<num_files;i++) {
            if(strcmp(label,strtok(line_array[i],"\t")) == 0) {
                printf("\t%d", atoi(strtok(NULL,"\t")));
            } else {
                if(end != 1) {
                    fprintf(stderr,"There was a label mismatch. I'll finish this line and exit");
                    end = 1;
                }
                printf("\tNA");
            }
        }
        printf("\n");
        if(end == 1) break;
    }

    //Clean things up
    for(i=0;i<num_files;i++) {
        free(line_array[i]);
        fclose(fp_array[i]);
    }
    free(line_array);
    free(fp_array);
    return 0;
}
