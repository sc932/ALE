// (C) 2010 Scott Clark, JGI, LBNL

// Compiling instructions
// $ cc -g -O2 ALEv2.c -o ALE -lz -lm
// $ ./ALE <parameter file>
// requires: ALE.h, kseq.h, zlib.h, geneTree.h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include "ALE.h" // included by geneTree.h
#include "geneTree.h"

struct SAMalignment{ // values in order of read in
    char readName[256];
    int outInfo;
    char refName[256];
    int mapStart;
    int mapPair;
    char cigar[256];
    char flag2[10];
    int mapEnd;
    int mapLen;
    char readSeq[256];
    char readQual[256];
    char XA[256];
    char MD[256];
    char NM[256];
    struct SAMalignment *pair;
};

typedef struct SAMalignment SAM_t;

int printSAM(SAM_t read){
    printf("SAM alignment:\n");
    printf("%s %i %s %i %i %s %s %i %i %s %s %s %s %s\n", read.readName, read.outInfo, read.refName, read.mapStart, read.mapPair, read.cigar, read.flag2, read.mapEnd, read.mapLen, read.readSeq, read.readQual, read.XA, read.MD, read.NM);
    
}

int hackedIntCast(char c){
    if(c == '0'){return 0;}
    if(c == '1'){return 1;}
    if(c == '2'){return 2;}
    if(c == '3'){return 3;}
    if(c == '4'){return 4;}
    if(c == '5'){return 5;}
    if(c == '6'){return 6;}
    if(c == '7'){return 7;}
    if(c == '8'){return 8;}
    if(c == '9'){return 9;}
}

double GetInsertProbNormal(const double sigma, const double point){
  //printf("Point: %f, p1: %f, p2: %f\n", point, erf((point + 0.5)/sqrt(2*sigma*sigma)), erf((point - 0.5)/sqrt(2*sigma*sigma)));
  return 0.5*(erf((abs(point) + 0.5)/sqrt(2*sigma*sigma)) - erf((abs(point) - 0.5)/sqrt(2*sigma*sigma)));
}

double getInsertLikelihood(SAM_t *read, float mu, float var){
    return GetInsertProbNormal(abs((float)read->mapLen) - mu, var);
}

// finds the likelihood of a string of misses in read from seqPos to matchLen
double likeMiss(SAM_t *read, int seqPos, int missLen){
    int i;
    double likelihood = 1.0;
    for(i = seqPos; i < seqPos + missLen; i++){
        likelihood = likelihood*((1.0 - QtoP[read->readQual[i] - 64])/3.0);
    }
    return likelihood;
}

// finds the likelihood of a string of matches in read from seqPos to matchLen
double likeMatch(SAM_t *read, int seqPos, int matchLen){
    int i;
    double likelihood = 1.0;
    for(i = seqPos; i < seqPos + matchLen; i++){
        likelihood = likelihood*(QtoP[read->readQual[i] - 64]);
    }
    return likelihood;
}

// takes in a read and returns the match likelihood (due to matches, mismatches, indels)
double getMatchLikelihood(SAM_t *read){
    
    int stop = 0;
    int pos = 5;
    int seqPos = 0;
    int base = 1;
    int matchLen, missLen, delLen, insLen;
    int totalMatch = 0, totalMiss = 0, totalDel = 0, totalIns = 0;
    int hardClipped = 0;
    double likelihood = 1.0;
    // parse MD field
    while(stop == 0){
        // matches
        matchLen = 0;
        while(isdigit(read->MD[pos])){
            matchLen = matchLen*10 + hackedIntCast(read->MD[pos]);
            pos++;
        }
        totalMatch += matchLen;
        likelihood = likelihood*likeMatch(read, seqPos, matchLen);
        seqPos += matchLen;
        // misses
        missLen = 0;
        while(read->MD[pos] == 'A' || read->MD[pos] == 'T' || read->MD[pos] == 'C' || read->MD[pos] == 'G' || read->MD[pos] == '0'){
            if(read->MD[pos] == 'A' || read->MD[pos] == 'T' || read->MD[pos] == 'C' || read->MD[pos] == 'G'){
                missLen++;
            }
            pos++;
        }
        totalMiss += missLen;
        likelihood = likelihood*likeMiss(read, seqPos, missLen);
        seqPos += missLen;
        // deletions
        delLen = 0;
        if(read->MD[pos] == '^'){
            pos++;
            while(isalpha(read->MD[pos])){
                delLen++;
                pos++;
            }
        }
        totalDel += delLen;
        // sees if we are at the end
        if(read->MD[pos] == '\0'){
            stop = 1;
        }
    }
    // parse CIGAR
    stop = 0;
    pos = 0;
    int num;
    while(stop == 0){
        // matches
        num = 0;
        while(isdigit(read->cigar[pos])){
            num = num*10 + hackedIntCast(read->cigar[pos]);
            pos++;
        }
        // insertions
        if(read->cigar[pos] == 'I'){
            totalIns += num;
        }
        pos++;
        // read to the end
        if(read->cigar[pos] == '\0'){
            stop = 1;
        }
    }
    if(totalDel + totalIns > 0){
        likelihood = 0.0; // can modify later to include insertions and deletions
    }
    printf("Found %i match(es).\n", totalMatch);
    printf("Found %i miss(es).\n", totalMiss);
    printf("Found %i deletion(s).\n", totalDel);
    printf("Found %i insertion(s).\n", totalIns);
    printf("Likelihood: %f.\n", likelihood);
    return likelihood;
}

int main(int argc, char **argv){
    // Handles input variations
    if (argc < 3) {
        if(argv[1][0] == '-' && argv[1][1] == 'h'){
            printf("%s", WELCOME_MSG);
            printf("Usage: %s [-options] input output\n\n", argv[0]);
            printf("%s", LONG_OPTIONS);
            return 0;
        }
        printf("Usage: %s [-options] input output\n", argv[0]);
        printf("%s", SHORT_OPTIONS);
        return 1;
    }
    
    // default parameter values
    float insertMeanInward = 200.0;
    float insertVarInward = 10.0;
    float insertMeanOutward = 200.0;
    float insertVarOutward = 10.0;
    
    // input and output files
    printf("Input file: %s\n", argv[argc - 2]);
    printf("Output file: %s\n", argv[argc - 1]);
    
    // attempt to open the input file
    FILE *ins = fopen(argv[argc - 2], "r");
    if(ins == NULL){
        printf("Error! Could not open input file: %s\n", argv[argc - 2]);
    }
    
    SAM_t read, readMate;
    double likelihoodRead1, likelihoodRead2, likelihoodInsert;
    // read in the first part of the read
    while( fscanf( ins, "%255s%i%255s%i%i%255s%10s%i%i%255s%255s%255s%255s%255s", &read.readName, &read.outInfo, &read.refName, &read.mapStart, &read.mapPair, &read.cigar, &read.flag2, &read.mapEnd, &read.mapLen, &read.readSeq, &read.readQual, &read.XA, &read.MD, &read.NM) > 0){
        
        printSAM(read);
        
        if (read.flag2[0] == '='){ // read in the mate, if it maps
             fscanf( ins, "%255s%i%255s%i%i%255s%10s%i%i%255s%255s%255s%255s%255s", &readMate.readName, &readMate.outInfo, &readMate.refName, &readMate.mapStart, &readMate.mapPair, &readMate.cigar, &readMate.flag2, &readMate.mapEnd, &readMate.mapLen, &readMate.readSeq, &readMate.readQual, &readMate.XA, &readMate.MD, &readMate.NM);
             
             printSAM(readMate);
             
             // compute the statitsics
             likelihoodRead1 = getMatchLikelihood(&read);
             likelihoodRead2 = getMatchLikelihood(&readMate);
             likelihoodInsert = getInsertLikelihood(&read, insertMeanInward, insertVarInward);
        }
    }
}