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
#include "ALEv2.h"
#include "ALElike.h"

int main(int argc, char **argv){
    // Handles input variations
    if (argc < 2) {
        printf("Usage: %s [-options] map assembly output\n", argv[0]);
        printf("%s", SHORT_OPTIONS);
        return 0;
    }
    if (argc < 4) {
        if(argc <= 1 || (argv[1][0] == '-' && argv[1][1] == 'h')){
            printf("%s", WELCOME_MSG);
            printf("Usage: %s [-options] map assembly output\n\n", argv[0]);
            printf("%s", LONG_OPTIONS);
            return 0;
        }else{
            printf("Usage: %s [-options] map assembly output\n", argv[0]);
            printf("%s", SHORT_OPTIONS);
            return 0;
        }
    }
    
    int options;
    int numberAssemblyPieces = 0;
    int kmerLen = 4;
    float insertLength = -1.0;
    float insertStd = -1.0;
    
    if(argc > 5) { // look for command line options
        for(options = 1; options < argc - 4; options++){ // search over all options
            if(strcmp(argv[options], "-nap") == 0){
                numberAssemblyPieces = atoi(argv[options+1]);
                options++;
            }else if(strcmp(argv[options], "-kmer") == 0){
                kmerLen = atoi(argv[options+1]);
                if(kmerLen > 6){
                    printf("-kmer option of %i not in range [2,6], set to default [4].\n", kmerLen);
                    kmerLen = 4;
                }
                options++;
            }else if(strcmp(argv[options], "-inl") == 0){
                insertLength = atof(argv[options+1]);
                if(insertLength <= 0){
                    printf("-inl option of %f not in range (0,inf), will be calculated from input.\n", insertLength);
                    insertLength = -1.0;
                }
                options++;
            }else if(strcmp(argv[options], "-ins") == 0){
                insertStd = atof(argv[options+1]);
                if(insertStd <= 0){
                    printf("-ins option of %f not in range (0,inf), will be calculated from input.\n", insertStd);
                    insertStd = -1.0;
                }
                options++;
            }else{
                printf("Could not find option %s\n", argv[options]);
            }
        }
    }
    
    // input and output files
    printf("Map file: %s\n", argv[argc - 3]);
    printf("Assembly file: %s\n", argv[argc - 2]);
    printf("Output file: %s\n", argv[argc - 1]);
    
    // attempt to open the input file
    FILE *ins = fopen(argv[argc - 3], "r");
    if(ins == NULL){
        printf("Error! Could not open map file: %s\n", argv[argc - 3]);
    }
    
    // attempt to open the input file
    gzFile *assemblyFile = gzopen(argv[argc - 2], "r");
    kseq_t *Aseq;
    if(assemblyFile == NULL){
        printf("Error! Could not open assembly file: %s\n", argv[argc - 2]);
    }
    
    // calculate the insert mean/std if not given
    if(insertLength == -1.0 || insertStd == -1.0){
        printf("Insert length and std not given, will be calculated from input map.\n");
        
    }
    
    printf("Reading in assembly...\n");
    Aseq = kseq_init(assemblyFile);
    
    if(numberAssemblyPieces == 0){
        numberAssemblyPieces = findNumAssemPieces(Aseq);
    }
    printf("Looking for %i assembly parts.\n", numberAssemblyPieces);
    gzclose(assemblyFile);
    assemblyFile = gzopen(argv[argc - 2], "r");
    Aseq = kseq_init(assemblyFile);
    
    assemblyT *theAssembly = malloc(sizeof(int) + (numberAssemblyPieces + 1)*sizeof(contig_t));
    theAssembly->contigs = malloc((numberAssemblyPieces + 1)*sizeof(contig_t));
    theAssembly->numContigs = numberAssemblyPieces;
    readAssembly(Aseq, theAssembly);
    printf("Done reading in assembly.\n");
    gzclose(assemblyFile);
    //printAssembly(theAssembly);
    
    SAM_t read, readMate;
    double likelihoodRead1, likelihoodRead2, likelihoodInsert;
    
    // initialize
    alignSet_t alignments;
    strcpy(alignments.name, "-1");
    alignments.nextAlignment = NULL;
    alignSet_t *currentAlignment = &alignments;
    alignSet_t *head = currentAlignment;
    alignSet_t *extension;
    
    //printAssembly(theAssembly);
    assemblySanityCheck(theAssembly);
    
    printf("Reading in the map and computing statistics...\n");
    // read in the first part of the read
    int keepGoing = 1;
    
    // calculate the insert mean/std if not given
    int GCtot = 0;
    double lengthTotal = 0.0;
    double lengthStd = 0.0;
    int readCount = 0;
    if(insertLength == -1 || insertStd == -1){
        int mapLens[mapLens_MAX], i;
        for(i = 0; i < mapLens_MAX; i++){
            mapLens[i] = 0;
        }
        int GCmaps[GCmaps_MAX];
        for(i = 0; i < GCmaps_MAX; i++){
            GCmaps[i] = 0;
        }
        printf("Insert length and std not given, will be calculated from input map.\n");
        
        while(keepGoing > 0){
            keepGoing = fscanf( ins, "%255s%i%255s%i%i%255s%10s%i%i%255s%255s%255s", read.readName, &read.outInfo, read.refName, &read.mapStart, &read.mapPair, read.cigar, read.flag2, &read.mapEnd, &read.mapLen, read.readSeq, read.readQual, read.XA);
            
            if(keepGoing < 1){break;}
            
            if(read.cigar[0] != '*'){ // see if there is an actual alignment there
                keepGoing = fscanf( ins, "%255s%255s", read.MD, read.NM);
            }else{
                strcpy(read.MD, "MD:Z:0");
                strcpy(read.NM, "NM:i:0");
            }
            
            if (read.flag2[0] == '=' || read.flag2[0] == '*'){ // read in the mate, if it maps
                keepGoing = fscanf( ins, "%255s%i%255s%i%i%255s%10s%i%i%255s%255s%255s", readMate.readName, &readMate.outInfo, readMate.refName, &readMate.mapStart, &readMate.mapPair, readMate.cigar, readMate.flag2, &readMate.mapEnd, &readMate.mapLen, readMate.readSeq, readMate.readQual, readMate.XA);

                if(keepGoing < 1){break;}
                
                if(read.cigar[0] != '*'){ // see if there is an actual alignment there
                    keepGoing = fscanf( ins, "%255s%255s", readMate.MD, readMate.NM);
                }else{
                    strcpy(readMate.MD, "MD:Z:0");
                    strcpy(readMate.NM, "NM:i:0");
                }
                
                GCtot = getGCtotal(read.readSeq, getSeqLen(read.readSeq), readMate.readSeq, getSeqLen(readMate.readSeq));
                lengthTotal += (double)abs(read.mapLen);
                mapLens[abs(read.mapLen)] += 1;
                GCmaps[GCtot] += 1;
                readCount++;
            }
        }
        insertLength = lengthTotal/(double)readCount;
        printf("Found sample avg insert length to be %f from %i mapped reads\n", insertLength, readCount);
        
        for(i = 0; i < mapLens_MAX; i++){
            if(mapLens[i] > 0){
                lengthStd += mapLens[i]*((double)i - insertLength)*((double)i - insertLength);
                //printf("i : mapLens[i] :: %i : %i\n", i, mapLens[i]);
            }
        }
        insertStd = sqrt(lengthStd/(double)(readCount-1));
        printf("Found sample length std to be %f\n", insertStd);
    }
    
    fclose(ins);
    ins = fopen(argv[argc - 3], "r");
    if(ins == NULL){
        printf("Error! Could not open map file: %s\n", argv[argc - 3]);
    }
    keepGoing = 1;
    while(keepGoing > 0){
        keepGoing = fscanf( ins, "%255s%i%255s%i%i%255s%10s%i%i%255s%255s%255s", read.readName, &read.outInfo, read.refName, &read.mapStart, &read.mapPair, read.cigar, read.flag2, &read.mapEnd, &read.mapLen, read.readSeq, read.readQual, read.XA);
        
        if(keepGoing < 1){break;}
        
        if(read.cigar[0] != '*'){ // see if there is an actual alignment there
            keepGoing = fscanf( ins, "%255s%255s", read.MD, read.NM);
        }else{
            strcpy(read.MD, "MD:Z:0");
            strcpy(read.NM, "NM:i:0");
        }
        
//         printf("Checking...\n");
//         keepGoing = assemblySanityCheck(theAssembly);
//         if(keepGoing < 1){break;}
//         
//         printSAM(read); // sanity check
        
        if (read.flag2[0] == '=' || read.flag2[0] == '*'){ // read in the mate, if it maps
            keepGoing = fscanf( ins, "%255s%i%255s%i%i%255s%10s%i%i%255s%255s%255s", readMate.readName, &readMate.outInfo, readMate.refName, &readMate.mapStart, &readMate.mapPair, readMate.cigar, readMate.flag2, &readMate.mapEnd, &readMate.mapLen, readMate.readSeq, readMate.readQual, readMate.XA);

            if(readMate.cigar[0] != '*'){ // see if there is an actual alignment there
                keepGoing = fscanf( ins, "%255s%255s", readMate.MD, readMate.NM);
            }else{
                strcpy(readMate.MD, "MD:Z:0");
                strcpy(readMate.NM, "NM:i:0");
            }
            
//             //printf("Checking...\n");
//             keepGoing = assemblySanityCheck(theAssembly);
//             if(keepGoing < 1){break;}
//              
//             printSAM(readMate); // sanity check
             
            // compute the statitsics
            likelihoodRead1 = getMatchLikelihood(&read);
            likelihoodRead2 = getMatchLikelihood(&readMate);
            likelihoodInsert = getInsertLikelihood(&read, insertLength, insertStd);
//             printf("Likelihoods: %12f %12f %12f\n", likelihoodRead1, likelihoodRead2, likelihoodInsert);
            
            if(read.cigar[0] == '*'){
                //printf("No alignment.\n");
            }else if(strcmp(currentAlignment->name, "-1") == 0){ // first alignment
                //printf("First alignment.\n");
                // copy in all the info
                strcpy(currentAlignment->name, read.readName);
                strcpy(currentAlignment->mapName, read.refName);
                currentAlignment->likelihood = likelihoodRead1*likelihoodRead2*likelihoodInsert;
                if(read.mapLen > 0){
                    currentAlignment->start1 = read.mapStart;
                    currentAlignment->end1 = read.mapStart + getSeqLen(read.readSeq);
                }else{
                    currentAlignment->start1 = read.mapStart - getSeqLen(read.readSeq);
                    currentAlignment->end1 = read.mapStart;
                }
                if(readMate.mapLen > 0){
                    currentAlignment->start2 = readMate.mapStart;
                    currentAlignment->end2 = readMate.mapStart + getSeqLen(readMate.readSeq);
                }else{
                    currentAlignment->start2 = readMate.mapStart - getSeqLen(readMate.readSeq);
                    currentAlignment->end2 = readMate.mapStart;
                }
            }else if(strcmp(currentAlignment->name, read.readName) == 0){ // test to see if this is another alignment of the current set or a new one
                // extend the set of alignments
                extension = malloc(sizeof(alignSet_t));
                currentAlignment->nextAlignment = extension;
                // copy in all the info
                strcpy(extension->name, read.readName);
                strcpy(extension->mapName, read.refName);
                extension->nextAlignment = NULL;
                extension->likelihood = likelihoodRead1*likelihoodRead2*likelihoodInsert;
                if(read.mapLen > 0){
                    extension->start1 = read.mapStart;
                    extension->end1 = read.mapStart + getSeqLen(read.readSeq);
                }else{
                    extension->start1 = read.mapStart - getSeqLen(read.readSeq);
                    extension->end1 = read.mapStart;
                }
                if(readMate.mapLen > 0){
                    extension->start2 = readMate.mapStart;
                    extension->end2 = readMate.mapStart + getSeqLen(readMate.readSeq);
                }else{
                    extension->start2 = readMate.mapStart - getSeqLen(readMate.readSeq);
                    extension->end2 = readMate.mapStart;
                }
                currentAlignment = extension;
                //printf("Same alignment!\n");
            }else{ // new alignment
                //printf("New alignment!\n");
                // do the statistics on *head, that read is exausted
//                 printAlignments(head);
                //printf("test\n");
                applyPlacement(head, theAssembly);
                // refresh head and current alignment
                head = currentAlignment;
                strcpy(currentAlignment->name, read.readName);
                strcpy(currentAlignment->mapName, read.refName);
                currentAlignment->likelihood = likelihoodRead1*likelihoodRead2*likelihoodInsert;
                if(read.mapLen > 0){
                    currentAlignment->start1 = read.mapStart;
                    currentAlignment->end1 = read.mapStart + getSeqLen(read.readSeq);
                }else{
                    currentAlignment->start1 = read.mapStart - getSeqLen(read.readSeq);
                    currentAlignment->end1 = read.mapStart;
                }
                if(readMate.mapLen > 0){
                    currentAlignment->start2 = readMate.mapStart;
                    currentAlignment->end2 = readMate.mapStart + getSeqLen(readMate.readSeq);
                }else{
                    currentAlignment->start2 = readMate.mapStart - getSeqLen(readMate.readSeq);
                    currentAlignment->end2 = readMate.mapStart;
                }
                currentAlignment->nextAlignment = NULL;
            }
        }
    }
    
    fclose(ins);
    
    //printAssembly(theAssembly);
    assemblySanityCheck(theAssembly);
    
    // clean up the final alignment
    applyPlacement(head, theAssembly);
    //printAlignments(head);
    
    printf("Done reading in the map.\n");
    
    // compute statistics on assembly
    printf("Computing k-mer statistics...\n");
    computeKmerStats(theAssembly, kmerLen);
    printf("Done computing k-mer statistics.\n");
    
    printf("Computing depth statistics...\n");
    computeDepthStats(theAssembly);
    printf("Done computing depth statistics.\n");
    
    FILE *out = fopen(argv[argc - 1], "w");
    if(out == NULL){
        printf("Error! Could not open output file: %s\n", argv[argc - 1]);
    }
    writeToOutput(theAssembly, out);
    
    printf("Done computing statistics.\nOutput is in file: %s\n", argv[argc - 1]);
    
    //printAssembly(theAssembly);
    //free(theAssembly);
}
