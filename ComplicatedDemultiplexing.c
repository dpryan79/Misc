#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int getFileNumber(char *line1, char *line2, char **bc1, char **bc2, int n) {
    int i, j, dist1=0;

    for(i=0; i<n; i++) {
        j = 0;
        dist1 = 0;
        while(bc1[i][j] != '\n' && bc1[i][j]) {
            if(bc1[i][j] != line1[j]) {
                dist1++;
                if(dist1 >2) break;
            }
            j++;
        }
//printf("1\t%s %s distance %i %s %s\n", line1, line2, dist1, bc1[i], bc2[i]);
        if(dist1>2) continue;
        j = 0;
        if(bc2[i]) {
            dist1 = 0;
            while(bc2[i][j] != '\n' && bc2[i][j]) {
                if(bc2[i][j] != line2[j]) {
                    dist1++;
                    if(dist1 >2) break;
                }
                j++;
            }
        }
//printf("2\t%s %s distance %i %s %s\n", line1, line2, dist1, bc1[i], bc2[i]);
        if(dist1 <= 2) {
//printf("%s %s matches %i %s %s\n", line1, line2, i, bc1[i], bc2[i]);
            return i;
        }
    }
//printf("%s %s matches nothing\n", line1, line2);
    return -1;
}

void processReads(FILE *in1, FILE *in2, FILE *in3, FILE *in4, FILE **out1, FILE **out2, char **bc1, char **bc2, int n, FILE *unmapped1, FILE *unmapped2) {
    char line11[1024], line12[1024], line13[1024], line14[1024];
    char line21[1024], line22[1024], line23[1024], line24[1024];
    char line[1024];
    int oNumber, i;
//    int foo = 0;

    while(fgets(line11, 1024, in2)) {
        fgets(line12, 1024, in2);
        fgets(line13, 1024, in2);
        fgets(line14, 1024, in2);
        fgets(line21, 1024, in3);
        fgets(line22, 1024, in3);
        fgets(line23, 1024, in3);
        fgets(line24, 1024, in3);

        line12[8] = '\0';
        line22[8] = '\0';
        oNumber = getFileNumber(line12, line22, bc1, bc2, n);
//printf("%s to file %i\n", line11, oNumber);
        if(oNumber >= 0) {
            for(i=0; i<4; i++) {
                fgets(line, 1024, in1);
                fputs(line, out1[oNumber]);
                fgets(line, 1024, in4);
                fputs(line, out2[oNumber]);
            }
        } else {
            for(i=0; i<4; i++) {
                fgets(line, 1024, in1);
                fputs(line, unmapped1);
                fgets(line, 1024, in4);
                fputs(line, unmapped2);
            }
        } //Else deal with unassigned stuff
//foo++;
//if(foo > 1000000) break;
    }
}

int main(int argc, char *argv[]) {
    FILE *in1=NULL, *in2=NULL, *in3=NULL, *in4=NULL, *sampleFile, *unmapped1, *unmapped2;
    FILE **outFiles1=NULL, **outFiles2=NULL;
    char **barcodes1=NULL, **barcodes2;
    int nSamples=0, i;
    char CMD[1024], line[1024], *project, *lib, *sample, *bc1, *bc2;

    if(argc != 6) {
        printf("Usage: %s Undetermined_S0_R1_001.fastq.gz Undetermined_S0_I1_001.fastq.gz Undetermined_S0_I2_001.fastq.gz Undetermined_S0_R2_001.fastq.gz samples.txt\n", argv[0]);
        printf("samples.txt is a csv with columns barcode1, barcode2, project, library, and sample name.\n");
        return 1;
    }

    //Open the input files
    sprintf(CMD, "zcat %s", argv[1]);
    in1 = popen(CMD, "r");
    sprintf(CMD, "zcat %s", argv[2]);
    in2 = popen(CMD, "r");
    sprintf(CMD, "zcat %s", argv[3]);
    in3 = popen(CMD, "r");
    sprintf(CMD, "zcat %s", argv[4]);
    in4 = popen(CMD, "r");
    sampleFile = fopen(argv[5], "r");

    while(fgets(line, 1024, sampleFile)) {
        bc1 = strtok(line, ",\n");
        bc2 = strtok(NULL, ",\n");
        project = strtok(NULL, ",\n");
        lib = strtok(NULL, ",\n");
        sample = strtok(NULL, ",\n");
        if(!sample) {
            sample = lib;
            lib = project;
            project = bc2;
            bc2 = NULL;
        }
        //Add everything to the lists
        nSamples++;
        outFiles1 = realloc(outFiles1, nSamples * sizeof(FILE*));
        outFiles2 = realloc(outFiles2, nSamples * sizeof(FILE*));
        barcodes1 = realloc(barcodes1, nSamples * sizeof(char*));
        barcodes2 = realloc(barcodes2, nSamples * sizeof(char*));
        sprintf(CMD, "mkdir -p BProject_%s/Sample_%s", project, lib);
        system(CMD);
        sprintf(CMD, "gzip > BProject_%s/Sample_%s/%s_R1.fastq.gz", project, lib, sample);
//printf("%i %s\n", nSamples-1, CMD);
        outFiles1[nSamples-1] = popen(CMD, "w");
        sprintf(CMD, "gzip > BProject_%s/Sample_%s/%s_R2.fastq.gz", project, lib, sample);
        outFiles2[nSamples-1] = popen(CMD, "w");
        barcodes1[nSamples-1] = strdup(bc1);
        if(bc2) {
            barcodes2[nSamples-1] = strdup(bc2);
        } else {
            barcodes2[nSamples-1] = NULL;
        }
    }
    fclose(sampleFile);

    sprintf(CMD, "gzip > BUnmapped_R1.fastq.gz");
    unmapped1 = popen(CMD, "w");
    sprintf(CMD, "gzip > BUnmapped_R2.fastq.gz");
    unmapped2 = popen(CMD, "w");

    processReads(in1, in2, in3, in4, outFiles1, outFiles2, barcodes1, barcodes2, nSamples, unmapped1, unmapped2);

    //Close things up
    pclose(in1);
    pclose(in2);
    pclose(in3);
    pclose(in4);
    for(i=0; i<nSamples; i++) {
        pclose(outFiles1[i]);
        pclose(outFiles2[i]);
        if(barcodes1[i]) free(barcodes1[i]);
        if(barcodes2[i]) free(barcodes2[i]);
    }
    free(barcodes1);
    free(barcodes2);
    free(outFiles1);
    free(outFiles2);
    free(unmapped1);
    free(unmapped2);

    return 0;
}
