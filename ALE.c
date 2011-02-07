// (C) 2010 Scott Clark, JGI, LBNL

// Compiling instructions
// $ cc -g -O2 ALE.c -o ALE -lz -lm
// $ ./ALE <parameter file>
// requires: ALE.h, kseq.h, zlib.h, geneTree.h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include "ALE.h"
#include "geneTree.h"
//#include "../MetassemblerC/MetassemblerCserial.h"

pairedRead_t FindBestPlacementsByIndex(pairedRead_t theRead, treeBranch_t *pRoot, assembly_t theAssembly, libInsOrProbs_t libProbs, int kmer);

int main(int argc, char **argv){
    int l, i, j; // for loops 
    time_t t0, t1; // for timing info
    int assemblyIdx;
    int readIdx;
    //int NUM_ASSEMBLY_PARTS = 0;
    int numReadFiles = 1;
    if (argc < 3) {
        printf("Usage: %s <parameter file> <output directory>\n", argv[0]);
    printf("%s", SHORT_OPTIONS);
        return 1;
    }
    //fp = gzopen(argv[1], "r");
    
    printf("Opening parameter file: %s\n", argv[1]);
    
    FILE *inputFile;
    
    inputFile = fopen(argv[1], "r");
    if( inputFile == NULL ){
      printf("Error opening parameter file!\n");
      return 0;
    }
    
    // open output up the files
    FILE *fo1, *fo2;
    
    pairedRead_t *theReads;
    
    FILE *cOut;
    char outputFile[256];
    outputFile[0] = '\0';
    strcat(outputFile, argv[2]);
    strcat(outputFile, "/ALEoutput.txt");
    //printf("OutputFile : %s\n", outputFile);
    cOut = fopen(outputFile, "w");
    if(cOut == NULL){
      printf("Error! Could not make output file!\n");
      printf("Please make sure that %s is a valid directory.\n", argv[2]);
      return 0;
    }
    
    char placementFile[256];
    char depthFile[256];
    placementFile[0] = '\0';
    depthFile[0] = '\0';
    strcat(placementFile, argv[2]);
    strcat(depthFile, argv[2]);
    strcat(placementFile, "/placementsOutput.txt");
    strcat(depthFile, "/depthOutput.txt");
    
    fo1 = fopen(placementFile, "w");
    fo2 = fopen(depthFile, "w");
    
    // read in the parameter file
    char bufferedInput[256];
    char tempInfo[64];
    char assemblyFile[256];
    unsigned char assemblyFileLen;
    int pos = 0;
    i = 0;
    char readFile[128][256];
    char stop = 0;
    int *readsInFile;
    int tempPointer;
    libInsOrProbs_t *libProbs;
    while(stop == 0){
      if(fgets(bufferedInput, 256, inputFile) == NULL){
    stop = 1;
      }else{
    //printf("%s", bufferedInput);
      if((char)bufferedInput[0] != '#'){ //comments
    //printf("Not on a comment line.\n");
    if(bufferedInput[0] == '>'){ // the assembly file
      //return 0;
      l = 1;
      while(bufferedInput[l] != '\n'){
        assemblyFile[l-1] = bufferedInput[l];
        l++;
      }
      assemblyFile[l-1] = '\0';
      assemblyFileLen = l;
      if(fgets(bufferedInput, 256, inputFile)){
        //printf("%s\n", bufferedInput);
      }else{printf("Error!\n");}
      numReadFiles = atoi(bufferedInput);
      //return 0;
      //readFile = malloc(numReadFiles*256*sizeof(char));
      readsInFile = malloc(numReadFiles*sizeof(int));
      libProbs = malloc(numReadFiles*sizeof(libInsOrProbs_t));
      
      //char readFile[numReadFiles][256];
      //int readsInFile[numReadFiles];
      //libInsOrProbs_t libProbs[numReadFiles];
    }else{ // read info
      l = 0;
      while(bufferedInput[l] != '\n'){
        readFile[i][l] = bufferedInput[l];
        l++;
      }
      readFile[i][l] = '\0';
      
      //printf("Reads file %i : %s\n", i, readFile[i]);
      
      if(fgets(bufferedInput, 256, inputFile)){
        //printf("%s\n", bufferedInput);
      }else{
        printf("Error!\n");
      } // read in the distribution info
      
      l = 0;
      pos = 0;
      //return 0;
      for(j = 0; j < 256; j++){
        //printf("char check: %c\n", bufferedInput[j]);
        if(bufferedInput[j] == '!'){
          j = 300;
        }
        while(bufferedInput[j] != ',' && bufferedInput[j] != '!'){
          tempInfo[l] = bufferedInput[j];
          l++;
          j++;
        }
        if(bufferedInput[j] == '!'){
          j = 300;
        }
        tempInfo[l] = '\0';
        l = 0;
        //printf("tempInfo holds: %s.\n", tempInfo);
        if(pos == 0){
          readsInFile[i] = atoi(tempInfo);
        }else if(pos == 1){
          libProbs[i].IOs[0].prob = atof(tempInfo);
        }else if(pos == 2){
          libProbs[i].IOs[0].mu = atof(tempInfo);
        }else if(pos == 3){
          libProbs[i].IOs[0].sigma = atof(tempInfo);
        }else if(pos == 4){
          libProbs[i].IOs[1].prob = atof(tempInfo);
        }else if(pos == 5){
          libProbs[i].IOs[1].mu = atof(tempInfo);
        }else if(pos == 6){
          libProbs[i].IOs[1].sigma = atof(tempInfo);
        }
        pos++;
      }
      i++;
    }
      }}
    }
    
    for(i = 0; i < numReadFiles; i++){
      printf("readsInFile: %i : %i\n", i, readsInFile[i]);
      printf("IOs[0]: %f, %f, %f\n", libProbs[i].IOs[0].prob, libProbs[i].IOs[0].mu, libProbs[i].IOs[0].sigma);
      printf("IOs[1]: %f, %f, %f\n", libProbs[i].IOs[1].prob, libProbs[i].IOs[1].mu, libProbs[i].IOs[1].sigma);
    }
    
    //return 0;
    
     
    
    /***************************************************************************
    ***** READ IN THE POTENTIAL ASSEMBLY/CONTIGS/SCAFFOLD **********************
    ***************************************************************************/
    
    // open up the assembly
    printf("Attemping to open assembly file: %s\n", assemblyFile);
    
    gzFile fp;
    kseq_t *seq;
    //return 0;
    //assemblyFile = "genData/genAssembly25.fna";
    fp = gzopen(&assemblyFile[0], "r");
    if( fp == NULL ){
      printf("Error opening input file!\n");
      return 0;
    }
    seq = kseq_init(fp);
    printf("Read in assembly file...\n");
    
    assembly_t theAssembly;
    i = 0;
    int NUM_ASSEMBLY_PARTS = 0;
    while ((l = kseq_read(seq)) >= 0) {
      
      NUM_ASSEMBLY_PARTS++;
      printf("NUM_ASSEMBLY_PARTS: %i\n", NUM_ASSEMBLY_PARTS);
      IncreaseAssemblyPartsByOne(&theAssembly, NUM_ASSEMBLY_PARTS);
      
      //l = kseq_read(seq); // read it in
      
      theAssembly.assemblyParts[i].seqLen = (int)(seq->seq.l); // set the length
      
      // set up the sequence memory
      theAssembly.assemblyParts[i].sequence = malloc(sizeof(char) * theAssembly.assemblyParts[i].seqLen);
      if(theAssembly.assemblyParts[i].sequence == NULL){
    printf("Memory error! (assembly sequence %i)\n", i);
    return 0;
      }
      
      // set up the depth calculations
      theAssembly.assemblyParts[i].depthInward = malloc(sizeof(float) * theAssembly.assemblyParts[i].seqLen);
      theAssembly.assemblyParts[i].depthOutward = malloc(sizeof(float) * theAssembly.assemblyParts[i].seqLen);
      if(theAssembly.assemblyParts[i].depthOutward == NULL){
    printf("Memory error! (assembly depth %i)\n", i);
    return 0;
      }
      
      for(j = 0; j < theAssembly.assemblyParts[i].seqLen; j++){
    theAssembly.assemblyParts[i].depthInward[0] = 0.0;
    theAssembly.assemblyParts[i].depthOutward[0] = 0.0;
      }
      
      // fill the sequence
      makeAssemblySeq(seq->seq.s, theAssembly.assemblyParts[i].seqLen, theAssembly.assemblyParts[i].sequence); // set the sequence
      i++;
    }
    
    
    // close up the file
    kseq_destroy(seq);
    //gzclose(fp); // for some reason causes a seg fault ???
    
//       printf("The assembly!\n");
//       for(i = 0; i < NUM_ASSEMBLY_PARTS; i++){
//         printf("Part %i:\n", i);
//         PrintAssembly(theAssembly.assemblyParts[i].sequence, theAssembly.assemblyParts[i].seqLen);
//       }
    
    
    printf("Building index tree...\n");
    t0 = time(NULL);
    
    int kmer = 10;
    treeBranch_t treeRoot, *pRoot = &treeRoot;
    treeRoot = MakeTree(theAssembly, kmer, NUM_ASSEMBLY_PARTS);
    printf("Done building tree. Took %i seconds.\n", (int)(time(NULL) - t0));
    t0 = time(NULL);
    
    /***************************************************************************
    **************** LOAD THE READS ********************************************
    ***************************************************************************/
    //return 0;
    double placementProb[NUM_ASSEMBLY_PARTS];
    //return 0;
    int numberPlacedOnPart[NUM_ASSEMBLY_PARTS];
    int unmappedReads = 0;
    double normalizerProb;
    float normalizer = 0.0;
    int startPos, endPos;
    int readFileIdx;
    
    int area1e, area1s, area2e, area2s;
    
    double depthAvgIn[NUM_ASSEMBLY_PARTS];
    double depthAvgOut[NUM_ASSEMBLY_PARTS];
    // log(prod(poisson)) = -lambda*n - log(lambda)*sum(k) - sum(log(k!))
    double poissonPart1in[NUM_ASSEMBLY_PARTS]; // sum(k)
    double poissonPart2in[NUM_ASSEMBLY_PARTS]; // sum(log(k!))
    double poissonPart1out[NUM_ASSEMBLY_PARTS]; // sum(k)
    double poissonPart2out[NUM_ASSEMBLY_PARTS]; // sum(log(k!))
    double logProbIn[NUM_ASSEMBLY_PARTS];
    double logProbOut[NUM_ASSEMBLY_PARTS];
    
    int qualityHistogram[63];
    int qualityHistogramNormalizer;
    
    double *insertLenHistogramIn;
    int insertLenHistogramInMax;
    double insertLenHistogramInNormalizer;
    double *insertLenHistogramOut;
    int insertLenHistogramOutMax;
    double insertLenHistogramOutNormalizer;
    
    double insertLenSampleMuIn[numReadFiles+1];
    double insertLenSampleSigmaIn[numReadFiles+1];
    double insertLenSampleMuOut[numReadFiles+1];
    double insertLenSampleSigmaOut[numReadFiles+1];
    
    fprintf(fo1, "Read Number, Placement number, total number of placements, likelihood, area1 start, area1 end, area 2 start, area 2 end, map type (0,1,4,5,8,9,12,13), assembly part\n");
    
    //return 0;
    int NUM_PAIRED_READS_ON_NODE;
    // START FOR LOOP OVER ALL THE READS FILES
    for(readFileIdx = 0; readFileIdx < numReadFiles; readFileIdx++){
      
      insertLenHistogramInMax = (int)(libProbs[readFileIdx].IOs[0].mu + 20*libProbs[readFileIdx].IOs[0].sigma);
      insertLenHistogramIn = malloc(sizeof(double) * insertLenHistogramInMax); // allow for up to a 20 sigma event
      insertLenHistogramOutMax = (int)(libProbs[readFileIdx].IOs[1].mu + 20*libProbs[readFileIdx].IOs[0].sigma);
      insertLenHistogramOut = malloc(sizeof(double) * insertLenHistogramOutMax); // allow for up to a 20 sigma event
      
      for(i = 0; i < 64; i++){
        qualityHistogram[i] = 0;
      }
      qualityHistogramNormalizer = 0;
      
      fprintf(fo1, "Library %i:\n", readFileIdx);
      fprintf(fo2, "Library %i:\n", readFileIdx);
      
      // Open up the reads
      printf("Attempting to open reads file: %s\n", readFile[readFileIdx]);
      fp = gzopen(readFile[readFileIdx], "r");
      seq = kseq_init(fp);
      printf("Read in reads file %s...\n", readFile[readFileIdx]);
      
      NUM_PAIRED_READS_ON_NODE = readsInFile[readFileIdx];
      printf("Attempting to read in %i reads.\n", NUM_PAIRED_READS_ON_NODE);
      // Read in the sequences to the paired_read structs
      //pairedRead_t theReads[NUM_PAIRED_READS_ON_NODE];
      theReads = malloc(sizeof(pairedRead_t)*NUM_PAIRED_READS_ON_NODE);
      i = 0;
      while ((l = kseq_read(seq)) >= 0) { // read in the sequence, if it exists
      //printf("name: %s\n", seq->name.s);

      theReads[i].readSubSeq[0].seqLen = (int)(seq->seq.l); // set the length
      // set the sequence
      makeSeq(seq->seq.s, theReads[i].readSubSeq[0].seqLen, theReads[i].readSubSeq[0].sequence);
      if (seq->qual.l){ // set the quality
        theReads[i].hasQval = 1;
        makeQual(seq->qual.s, theReads[i].readSubSeq[0].seqLen, theReads[i].readSubSeq[0].qval, qualityHistogram);
        qualityHistogramNormalizer += theReads[i].readSubSeq[0].seqLen;
      }else{
        makeQualPerfect(theReads[i].readSubSeq[0].seqLen, theReads[i].readSubSeq[0].qval);
        theReads[i].hasQval = 0;
      }
      
      l = kseq_read(seq); // read in the next sequence
      
      theReads[i].readSubSeq[1].seqLen = (int)(seq->seq.l); // set the length
      // set the sequence
      makeSeq(seq->seq.s, theReads[i].readSubSeq[1].seqLen, theReads[i].readSubSeq[1].sequence);
      if (seq->qual.l){ // set the quality
        theReads[i].hasQval = 1;
        makeQual(seq->qual.s, theReads[i].readSubSeq[1].seqLen, theReads[i].readSubSeq[1].qval, qualityHistogram);
        qualityHistogramNormalizer += theReads[i].readSubSeq[1].seqLen;
      }else{
        makeQualPerfect(theReads[i].readSubSeq[1].seqLen, theReads[i].readSubSeq[1].qval);
        theReads[i].hasQval = 0;
      }
      
      theReads[i].numPlacements = 0;
      i++; // go to the next read
      }
      // close up the file
      kseq_destroy(seq);
      gzclose(fp);
      
  /*******************************************************************************
  ************* PRINTING CHECK ***************************************************
  *******************************************************************************/

      // print out the sequences
      printf("The sequences!\n");
      char seqer[4];
      for(i = 0; i < NUM_PAIRED_READS_ON_NODE; i++){
          printf("Sequence %i:\n",i);
          printf("Read 1(%i):\n",i);
          PrintSequence(theReads[i].readSubSeq[0].sequence, theReads[i].readSubSeq[0].seqLen);
  	PrintQuality(theReads[i].readSubSeq[0].qval, theReads[i].readSubSeq[0].seqLen);
  	printf("Read 2(%i):\n",i);
          PrintSequence(theReads[i].readSubSeq[1].sequence, theReads[i].readSubSeq[1].seqLen);
  	PrintQuality(theReads[i].readSubSeq[1].qval, theReads[i].readSubSeq[1].seqLen);
  	//printf("%f", getQuality(theReads[i].readSubSeq[1].qval, 4));
      }
      

      printf("Done reading. Took %i seconds.\n", (int)(time(NULL) - t0));
      
//       // find average read score (see writeup/paper)
//       
//       double readLenMu = 85.0;
//       double readLenSigma = 7.0;
//       double averageReadScore = 0.0;
//       double averagePoissonScore[NUM_ASSEMBLY_PARTS];
//       double arsTlikelihood = 0.0;
//       double pQ, pL, pR1, pR2, pO;
//       int arsQ, arsInsertLen, arsRead1len, arsRead2len, arsOrient, arsDepth;
//       double arsThresh = 0.000000000001;
//       for(arsQ = 0; arsQ < 64; arsQ++){ // quality loop
// 	pQ = (float)qualityHistogram[arsQ]/(float)qualityHistogramNormalizer; // prob of this quality
// 	if(pQ > arsThresh){ // threshold test
// 	  for(arsOrient = 0; arsOrient < 2; arsOrient++){
// 	    pO = libProbs[readFileIdx].IOs[arsOrient].prob; // prob of this orientation
// 	    if(pQ*pO > arsThresh){ // threshold test
// 	      for(arsInsertLen = 0; arsInsertLen < 2*libProbs[readFileIdx].IOs[arsOrient].mu; arsInsertLen++){
// 		pL = GetInsertProbNormal(libProbs[readFileIdx].IOs[arsOrient].sigma, (float)arsInsertLen - libProbs[readFileIdx].IOs[arsOrient].mu);
// 		if(pQ*pO*pL > arsThresh){ // threshold test
// 		  for(arsRead1len = 10; arsRead1len < arsInsertLen; arsRead1len++){
// 		    pR1 = GetInsertProbNormal(readLenSigma, arsRead1len - readLenMu);
// 		    if(pQ*pO*pL*pR1 > arsThresh){ // threshold test
// 		      for(arsRead2len = 10; arsRead2len < arsInsertLen; arsRead2len++){
// 			pR2 = GetInsertProbNormal(readLenSigma, arsRead2len - readLenMu);
// 			if(pQ*pO*pL*pR1*pR2 > arsThresh){ // threshold test
// 			  arsTlikelihood = log(QtoP[arsQ]*QtoP[arsQ] + (1-QtoP[arsQ])*(1-QtoP[arsQ])/3.0)*(float)(arsRead1len + arsRead2len) + log(pL) + log(pO);
// 			  averageReadScore += pQ*pO*pL*pR1*pR2*arsTlikelihood;
// 			}
// 		      }
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }
//       printf("AverageScoreLikelihood: %f\n", averageReadScore);
      
      
      printf("Finding placements...\n");
      
      /***************************************************************************
      ************** FIND PLACEMENTS *********************************************
      ***************************************************************************/
      
      //placement_t tempPlacements[N_PLACEMENTS];
      //int numTempPlacements = 0;
      
      for(readIdx = 0; readIdx < NUM_PAIRED_READS_ON_NODE; readIdx++){ // cycle over all the reads
    
    // print out timing info ever 100,000 reads
        if(readIdx%1 == 0){
          t1 = time(NULL);
          printf("Time elapsed: %f\n", difftime(t1,t0));
          t0 = t1;
          printf("Checking read %i\n", readIdx);
        }
        // find the best placements for the reads using an prefix of length <kmer> index search first
        theReads[readIdx] = FindBestPlacementsByIndex(theReads[readIdx], pRoot, theAssembly, libProbs[readFileIdx], kmer);
        printf("Checked.\n");
      }
      printf("Time elapsed: %f\n", difftime(t1,t0));
      t0 = t1;
      
      /***************************************************************************
      ************** STATS PART **************************************************
      ***************************************************************************/
      
      
      
      
      unmappedReads = 0;
      for(i = 0; i < NUM_ASSEMBLY_PARTS; i++){
    placementProb[i] = 0.0;
    numberPlacedOnPart[i] = 0;
      }
      
      for(readIdx = 0; readIdx < NUM_PAIRED_READS_ON_NODE; readIdx++){
    //fprintf(fo1, "Read %i, printing %i placements:\n", readIdx, theReads[readIdx].numPlacements);
    normalizerProb = 0.0;
    if(theReads[readIdx].numPlacements == 0){
      unmappedReads++;
    }
    for(i = 0; i < theReads[readIdx].numPlacements; i++){
      //fprintf(fo1, "L: %f, o1: %i, o2: %i, i: %i an: %i\n", theReads[readIdx].placements[i].likelihood, theReads[readIdx].placements[i].offset1, theReads[readIdx].placements[i].offset2, theReads[readIdx].placements[i].placeInfo, theReads[readIdx].placements[i].assemPart);
      area1s = 0;
      area1e = 0;
      area2s = 0;
      area2e = 0;
      if(theReads[readIdx].placements[i].placeInfo == 13 || theReads[readIdx].placements[i].placeInfo == 9){
        area1s = theReads[readIdx].placements[i].offset1;
        area1e = theReads[readIdx].placements[i].offset1 + theReads[readIdx].readSubSeq[0].seqLen;
        area2s = theReads[readIdx].placements[i].offset2 + kmer - 1 - theReads[readIdx].readSubSeq[1].seqLen;
        area2e = theReads[readIdx].placements[i].offset2 + kmer - 1;
      }else if(theReads[readIdx].placements[i].placeInfo == 12 || theReads[readIdx].placements[i].placeInfo == 8){
        area1s = theReads[readIdx].placements[i].offset2;
        area1e = theReads[readIdx].placements[i].offset2 + theReads[readIdx].readSubSeq[1].seqLen;
        area2s = theReads[readIdx].placements[i].offset1 + kmer - 1 - theReads[readIdx].readSubSeq[0].seqLen;
        area2e = theReads[readIdx].placements[i].offset1 + kmer - 1;
      }else if(theReads[readIdx].placements[i].placeInfo == 5 || theReads[readIdx].placements[i].placeInfo == 1){
        area1s = theReads[readIdx].placements[i].offset1 + kmer - 1 - theReads[readIdx].readSubSeq[0].seqLen;
        area1e = theReads[readIdx].placements[i].offset1 + kmer - 1;
        area2s = theReads[readIdx].placements[i].offset2;
        area2e = theReads[readIdx].placements[i].offset2 + theReads[readIdx].readSubSeq[1].seqLen;
      }else{
        area1s = theReads[readIdx].placements[i].offset2 + kmer - 1 - theReads[readIdx].readSubSeq[1].seqLen;
        area1e = theReads[readIdx].placements[i].offset2 + kmer - 1;
        area2s = theReads[readIdx].placements[i].offset1;
        area2e = theReads[readIdx].placements[i].offset1 + theReads[readIdx].readSubSeq[0].seqLen;
      }
      fprintf(fo1, "%i,%i,%i,%.12f,%i,%i,%i,%i,%i,%i\n", readIdx, i, theReads[readIdx].numPlacements, theReads[readIdx].placements[i].likelihood, area1s, area1e, area2s, area2e, theReads[readIdx].placements[i].placeInfo, theReads[readIdx].placements[i].assemPart);
      normalizerProb += theReads[readIdx].placements[i].likelihood;
    }
    theReads[readIdx].placementNormalizer = normalizerProb;
    for(i = 0; i < theReads[readIdx].numPlacements; i++){
      if(i == 0){ // only place the best one
        if(theReads[readIdx].placements[i].likelihood > LIKELIHOOD_THRESHOLD){
          placementProb[theReads[readIdx].placements[i].assemPart] += log(theReads[readIdx].placements[i].likelihood);
          // theReads[readIdx].placements[i].likelihood/normalizerProb*
          //placementProb[theReads[readIdx].placements[i].assemPart] += log(theReads[readIdx].placements[i].likelihood);
          numberPlacedOnPart[theReads[readIdx].placements[i].assemPart] += 1;
        }
      }
    }
// 	printf("Read %i, printing %i placements:\n", readIdx, theReads[readIdx].numPlacements);
// 	PrintPlacements(theReads[readIdx]);
      }
      
//       printf("Time elapsed (Stats 1): %f\n", difftime(t1,t0));
//       t0 = t1;
//       
//       // find the depths
//       
//       double assemPartNormalizer[NUM_ASSEMBLY_PARTS];
//       
//       for(readIdx = 0; readIdx < NUM_PAIRED_READS_ON_NODE; readIdx++){
// 	normalizer = 0.0;
// 	for(i = 0; i < NUM_ASSEMBLY_PARTS; i++){
// 	  assemPartNormalizer[i] = 0.0;
// 	}
// 	for(i = 0; i < theReads[readIdx].numPlacements; i++){
// 	  normalizer += theReads[readIdx].placements[i].likelihood;
// 	  assemPartNormalizer[theReads[readIdx].placements[i].assemPart] += theReads[readIdx].placements[i].likelihood; // want to normalize accross different assembly parts
// 	}
// 	//return 0;
// 	for(i = 0; i < theReads[readIdx].numPlacements; i++){
// 	  if((theReads[readIdx].placements[i].placeInfo == 9) || (theReads[readIdx].placements[i].placeInfo == 13)){ // oriented inward, left sequence first
// 	    //L part
// 	    startPos = intMax(0, theReads[readIdx].placements[i].offset1);
// 	    endPos = intMin(startPos + theReads[readIdx].readSubSeq[0].seqLen, theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].seqLen);
// 	    for(j = startPos; j < endPos; j++){ // first part
// 	      theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].depthInward[j] += (float)theReads[readIdx].placements[i].likelihood/normalizer*assemPartNormalizer[theReads[readIdx].placements[i].assemPart]/normalizer;
// 	    }
// 	    // R part
// 	    endPos = intMin(startPos + theReads[readIdx].readSubSeq[1].seqLen, theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].seqLen);
// 	    // insert len histogram change:
// 	    if(endPos - startPos < insertLenHistogramInMax){
// 	      insertLenHistogramIn[endPos - startPos] += theReads[readIdx].placements[i].likelihood/theReads[readIdx].placementNormalizer;
// 	    }
// 	    startPos = intMax(0, theReads[readIdx].placements[i].offset2 + kmer - theReads[readIdx].readSubSeq[1].seqLen);
// 	    for(j = startPos; j < endPos; j++){ // second part
// 	      theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].depthInward[j] += (float)theReads[readIdx].placements[i].likelihood/normalizer*assemPartNormalizer[theReads[readIdx].placements[i].assemPart]/normalizer;
// 	    }
// 	  }else if((theReads[readIdx].placements[i].placeInfo == 0) || (theReads[readIdx].placements[i].placeInfo == 4)){ // oriented outward, Right seq first
// 	    //R part
// 	    startPos = intMax(0, theReads[readIdx].placements[i].offset2 + kmer - theReads[readIdx].readSubSeq[1].seqLen);
// 	    endPos = intMin(startPos + theReads[readIdx].readSubSeq[0].seqLen, theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].seqLen);
// 	    for(j = startPos; j < endPos; j++){ // first part
// 	      theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].depthOutward[j] += (float)theReads[readIdx].placements[i].likelihood/normalizer*assemPartNormalizer[theReads[readIdx].placements[i].assemPart]/normalizer;
// 	    }
// 	    // L part
// 	    endPos = intMin(startPos + theReads[readIdx].readSubSeq[1].seqLen ,theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].seqLen);
// 	    // insert len histogram change:
// 	    if(endPos - startPos < insertLenHistogramOutMax){
// 	      insertLenHistogramOut[endPos - startPos] += theReads[readIdx].placements[i].likelihood/theReads[readIdx].placementNormalizer;
// 	    }
// 	    startPos = intMax(0, theReads[readIdx].placements[i].offset1);
// 	    for(j = startPos; j < endPos; j++){ // second part
// 	      theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].depthOutward[j] += (float)theReads[readIdx].placements[i].likelihood/normalizer*assemPartNormalizer[theReads[readIdx].placements[i].assemPart]/normalizer;
// 	    }
// 	  }else if((theReads[readIdx].placements[i].placeInfo == 1) || (theReads[readIdx].placements[i].placeInfo == 5)){// oriented inward, right sequence first
// 	    //R part
// 	    startPos = intMax(0, theReads[readIdx].placements[i].offset2);
// 	    endPos = intMin(theReads[readIdx].placements[i].offset2 + theReads[readIdx].readSubSeq[1].seqLen, theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].seqLen);
// 	    for(j = startPos; j < endPos; j++){ // first part
// 	      theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].depthInward[j] += (float)theReads[readIdx].placements[i].likelihood/normalizer*assemPartNormalizer[theReads[readIdx].placements[i].assemPart]/normalizer;
// 	    }
// 	    // L part
// 	    endPos = intMin(startPos + theReads[readIdx].readSubSeq[0].seqLen, theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].seqLen);
// 	    // insert len histogram change:
// 	    if(endPos - startPos < insertLenHistogramInMax){
// 	      insertLenHistogramIn[endPos - startPos] += theReads[readIdx].placements[i].likelihood/theReads[readIdx].placementNormalizer;
// 	    }
// 	    startPos = intMax(0, theReads[readIdx].placements[i].offset1 + kmer - theReads[readIdx].readSubSeq[0].seqLen);
// 	    for(j = startPos; j < endPos; j++){ // second part
// 	      theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].depthInward[j] += (float)theReads[readIdx].placements[i].likelihood/normalizer*assemPartNormalizer[theReads[readIdx].placements[i].assemPart]/normalizer;
// 	    }
// 	  }else{// oriented outward, left sequence first
// 	    //L part
// 	    startPos = theReads[readIdx].placements[i].offset1 + kmer - theReads[readIdx].readSubSeq[1].seqLen;
// 	    endPos = intMin(startPos + theReads[readIdx].readSubSeq[1].seqLen, theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].seqLen);
// 	    for(j = startPos; j < endPos; j++){ // first part
// 	      theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].depthOutward[j] += (float)theReads[readIdx].placements[i].likelihood/normalizer*assemPartNormalizer[theReads[readIdx].placements[i].assemPart]/normalizer;
// 	    }
// 	    // R part
// 	    endPos = intMin(startPos + theReads[readIdx].readSubSeq[0].seqLen, theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].seqLen);
// 	    // insert len histogram change:
// 	    if(endPos - startPos < insertLenHistogramOutMax){
// 	      insertLenHistogramOut[endPos - startPos] += theReads[readIdx].placements[i].likelihood/theReads[readIdx].placementNormalizer;
// 	    }
// 	    startPos = theReads[readIdx].placements[i].offset2;
// 	    for(j = startPos; j < endPos; j++){ // second part
// 	      theAssembly.assemblyParts[theReads[readIdx].placements[i].assemPart].depthOutward[j] += (float)theReads[readIdx].placements[i].likelihood/normalizer*assemPartNormalizer[theReads[readIdx].placements[i].assemPart]/normalizer;
// 	    }
// 	  }
// 	}
//       }
//       
//       printf("Time elapsed (Stats 2): %f\n", difftime(t1,t0));
//       t0 = t1;
//       
//       //calculate the placement sample insert mean, mu
//       insertLenSampleMuIn[readFileIdx] = 0.0;
//       insertLenSampleSigmaIn[readFileIdx] = 0.0;
//       insertLenSampleMuOut[readFileIdx] = 0.0;
//       insertLenSampleSigmaOut[readFileIdx] = 0.0;
//       insertLenHistogramInNormalizer = 0.0;
//       insertLenHistogramOutNormalizer = 0.0;
//       //mean
//       for(i = 0; i < insertLenHistogramInMax-1; i++){
// 	if(insertLenHistogramIn[i] > 0.1){
// 	  //printf("i : hist = %i : %f\n", i, insertLenHistogramIn[i]);
// 	}
// 	insertLenHistogramInNormalizer += insertLenHistogramIn[i];
// 	insertLenSampleMuIn[readFileIdx] += (double)i*insertLenHistogramIn[i];
//       }
//       insertLenSampleMuIn[readFileIdx] = insertLenSampleMuIn[readFileIdx]/(double)insertLenHistogramInNormalizer;
//       for(i = 0; i < insertLenHistogramOutMax-1; i++){
// 	insertLenHistogramOutNormalizer += insertLenHistogramOut[i];
// 	insertLenSampleMuOut[readFileIdx] += (double)i*insertLenHistogramOut[i];
//       }
//       insertLenSampleMuOut[readFileIdx] = insertLenSampleMuOut[readFileIdx]/(double)insertLenHistogramOutNormalizer;
//       //std
//       for(i = 0; i < insertLenHistogramInMax-1; i++){
// 	insertLenSampleSigmaIn[readFileIdx] += ((double)i-insertLenSampleMuIn[readFileIdx])*((double)i-insertLenSampleMuIn[readFileIdx])*insertLenHistogramIn[i];
//       }
//       insertLenSampleSigmaIn[readFileIdx] = sqrt(insertLenSampleSigmaIn[readFileIdx]/(double)insertLenHistogramInNormalizer);
//       for(i = 0; i < insertLenHistogramOutMax-1; i++){
// 	insertLenSampleSigmaOut[readFileIdx] += ((double)i-insertLenSampleMuOut[readFileIdx])*((double)i-insertLenSampleMuOut[readFileIdx])*insertLenHistogramOut[i];
//       }
//       insertLenSampleSigmaOut[readFileIdx] = sqrt(insertLenSampleSigmaOut[readFileIdx]/(double)insertLenHistogramOutNormalizer);
//       
//       printf("Sample inward mu, sigma: %f, %f (%f samples)\n", insertLenSampleMuIn[readFileIdx], insertLenSampleSigmaIn[readFileIdx], insertLenHistogramInNormalizer);
//       printf("Sample outward mu, sigma: %f, %f (%f samples)\n", insertLenSampleMuOut[readFileIdx], insertLenSampleSigmaOut[readFileIdx], insertLenHistogramOutNormalizer);
//       
//       free(insertLenHistogramIn);
//       free(insertLenHistogramOut);
//       
//       printf("Time elapsed (Stats 3): %f\n", difftime(t1,t0));
//       t0 = t1;
//       
//       // print out the depths
//       for(j = 0; j < NUM_ASSEMBLY_PARTS; j++){
// 	depthAvgIn[j] = 0.0;
// 	depthAvgOut[j] = 0.0;
// 	poissonPart1in[j] = 0.0;
// 	poissonPart2in[j] = 0.0;
// 	poissonPart1out[j] = 0.0;
// 	poissonPart2out[j] = 0.0;
// 	for(i = 0; i < theAssembly.assemblyParts[j].seqLen; i+=POISSON_COMB){
// 	  fprintf(fo2, "%i,%i,%f,%f\n", j, i, theAssembly.assemblyParts[j].depthInward[i], theAssembly.assemblyParts[j].depthOutward[i]);
// 	  if((i > libProbs[readFileIdx].IOs[0].mu) && (i < theAssembly.assemblyParts[j].seqLen - libProbs[readFileIdx].IOs[0].mu)){ // cut off the ends
// 	    depthAvgIn[j] += theAssembly.assemblyParts[j].depthInward[i]/(float)floor(theAssembly.assemblyParts[j].seqLen - 2*libProbs[readFileIdx].IOs[0].mu);
// 	    poissonPart1in[j] += theAssembly.assemblyParts[j].depthInward[i];
// 	    //poissonPart2in[j] += lnfact((double)theAssembly.assemblyParts[j].depthInward[i] - 1);
// 	    poissonPart2in[j] += lgamma((double)theAssembly.assemblyParts[j].depthInward[i] + 1);
// 	  }
// 	  if((i > libProbs[readFileIdx].IOs[1].mu) && (i < theAssembly.assemblyParts[j].seqLen - libProbs[readFileIdx].IOs[1].mu)){
// 	    depthAvgOut[j] += theAssembly.assemblyParts[j].depthOutward[i]/(float)floor(theAssembly.assemblyParts[j].seqLen - 2*libProbs[readFileIdx].IOs[1].mu);
// 	    poissonPart1out[j] += theAssembly.assemblyParts[j].depthOutward[i];
// 	    //poissonPart2out[j] += lnfact((double)theAssembly.assemblyParts[j].depthOutward[i] - 1);
// 	    poissonPart2out[j] += lgamma((double)theAssembly.assemblyParts[j].depthOutward[i] + 1);
// 	  }
// 	//printf("%i: %f + %f = %f\n", i, theAssembly.assemblyParts[0].depthInward[i], theAssembly.assemblyParts[0].depthOutward[i], theAssembly.assemblyParts[0].depthInward[i] + theAssembly.assemblyParts[0].depthOutward[i]);
// 	}
// 	
//   //       logProbIn[j] = 0;
//   //       for(i = 0; i < theAssembly.assemblyParts[j].seqLen; i++){
//   // 	if((i > libProbs.IOs[0].mu) && (i < theAssembly.assemblyParts[j].seqLen - libProbs.IOs[0].mu)){
//   // 	  logProbIn[j] += -depthAvgIn[j] + log(depthAvgIn[j])*theAssembly.assemblyParts[j].depthInward[i] - lgamma((double)theAssembly.assemblyParts[j].depthInward[i] + 1);
//   // 	  printf("log %i: %f\n", i, -depthAvgIn[j] + log(depthAvgIn[j])*theAssembly.assemblyParts[j].depthInward[i] - lgamma((double)theAssembly.assemblyParts[j].depthInward[i] + 1));
//   // 	}
//   //       }
//   //       printf("Poisson1 %i: %f, %f\n", j, poissonPart1in[j], poissonPart2in[j]);
//   //       printf("All the infos: %i, %f, %f, %f, %f, %f\n", j, depthAvgIn[j], (float)theAssembly.assemblyParts[j].seqLen - 2*libProbs.IOs[0].mu, log(depthAvgIn[j]), poissonPart1in[j], poissonPart2in[j]);
// 	// poisson pdf = lambda^k*e^-lambda / k!
// 	// log pp = -lambda + k log lambda - log k!
// 	averagePoissonScore[j] = 0.0;
// 	//printf("mu = %f\n", depthAvgIn[j]);
// 	for(arsDepth = 2; arsDepth < 2*depthAvgIn[j]; arsDepth++){
// 	  if(poissonInt(arsDepth, depthAvgIn[j]) > 0.00000001){
// 	    //printf("Tests %i: %f * %f = %f\n", arsDepth, log(poissonInt(arsDepth, depthAvgIn[j])), poissonInt(arsDepth, depthAvgIn[j]), log(poissonInt(arsDepth, depthAvgIn[j]))*poissonInt(arsDepth, depthAvgIn[j]));
// 	    
// 	    averagePoissonScore[j] += log(poissonInt(arsDepth, depthAvgIn[j]))*poissonInt(arsDepth, depthAvgIn[j]);
// 	  }
// 	}
// 	printf("AveragePoissonScore: %f\n", averagePoissonScore[j]);
// 	logProbIn[j] = -1.0*depthAvgIn[j]*(double)((float)theAssembly.assemblyParts[j].seqLen - 2*libProbs[readFileIdx].IOs[0].mu) + log(depthAvgIn[j])*poissonPart1in[j] - poissonPart2in[j];
// 	logProbOut[j] = -1.0*depthAvgOut[j]*(double)((float)theAssembly.assemblyParts[j].seqLen - 2*libProbs[readFileIdx].IOs[1].mu) + log(depthAvgOut[j])*poissonPart1out[j] - poissonPart2out[j];
// 	printf("Poisson %i: %f, %f, %f\n", j, logProbIn[j], logProbOut[j], (logProbIn[j] + logProbOut[j])/((float)theAssembly.assemblyParts[j].seqLen - 2*libProbs[readFileIdx].IOs[0].mu));
// 	printf("Placement %i: %f, %f\n", j, placementProb[j], placementProb[j]/(double)numberPlacedOnPart[j]);
// 	printf("With %i unmapped reads of %i: %f\n", unmappedReads, NUM_PAIRED_READS_ON_NODE, (double)unmappedReads/(double)NUM_PAIRED_READS_ON_NODE);
// 	fprintf(cOut, "Library %i:\n", readFileIdx);
// 	fprintf(cOut, "Poisson %i: %f, %f, %f\n", j, logProbIn[j], logProbOut[j], (logProbIn[j] + logProbOut[j])/((float)theAssembly.assemblyParts[j].seqLen - 2*libProbs[readFileIdx].IOs[0].mu));
// 	fprintf(cOut, "Placement %i: %f, %f\n", j, placementProb[j], placementProb[j]/(double)numberPlacedOnPart[j]);
// 	fprintf(cOut, "With %i unmapped reads of %i: %f\n", unmappedReads, NUM_PAIRED_READS_ON_NODE, (double)unmappedReads/(double)NUM_PAIRED_READS_ON_NODE);
//       }
//       
//       printf("Time elapsed (Stats 4): %f\n", difftime(t1,t0));
//       t0 = t1;

      free(theReads);
    
    }
    // END FOR LOOP OVER ALL THE READS FILES
    
    fclose(fo1);
    fclose(fo2);
    fclose(cOut);
    
    // causes a segmentation fault for some reason...
//     for(i = 0; i < NUM_ASSEMBLY_PARTS; i++){
//       //free(theAssembly.assemblyParts[i].sequence);
//       //free(theAssembly.assemblyParts[i].depthOutward);
//       //free(theAssembly.assemblyParts[i].depthInward);
//     }
    
}




pairedRead_t FindBestPlacementsByIndex(pairedRead_t theRead, treeBranch_t *pRoot, assembly_t theAssembly, libInsOrProbs_t libProbs, int kmer){
  // build the prefixes
  char prefix1[kmer], bprefix1[kmer]; // prefix forward and backward first part of read
  char prefix1c[kmer], bprefix1c[kmer]; // compliments
  char prefix2[kmer], bprefix2[kmer]; // prefix forward and backward first part of read
  char prefix2c[kmer], bprefix2c[kmer]; // compliments
  int i, j;

  for(i = 0; i < kmer; i++){
    prefix1[i] = getCharFromSeqByLoc(theRead.readSubSeq[0].sequence, i);
    prefix1c[i] = getComplimentRes(getCharFromSeqByLoc(theRead.readSubSeq[0].sequence, i));
    bprefix1[i] = getCharFromSeqByLoc(theRead.readSubSeq[0].sequence, kmer - 1 - i);
    bprefix1c[i] = getComplimentRes(getCharFromSeqByLoc(theRead.readSubSeq[0].sequence, kmer - 1 - i));
    prefix2[i] = getCharFromSeqByLoc(theRead.readSubSeq[1].sequence, i);
    prefix2c[i] = getComplimentRes(getCharFromSeqByLoc(theRead.readSubSeq[1].sequence, i));
    bprefix2[i] = getCharFromSeqByLoc(theRead.readSubSeq[1].sequence, kmer - 1 - i);
    bprefix2c[i] = getComplimentRes(getCharFromSeqByLoc(theRead.readSubSeq[1].sequence, kmer - 1 - i));
  }
  
//   //printf("prefix: %.10s\n", prefix1);
//   //printf("prefix1c: %.10s\n", prefix1c);
//   //printf("bprefix1: %.10s\n", bprefix1);
//   //printf("bprefix1c: %.10s\n", bprefix1c);
//   //printf("prefix2: %.10s\n", prefix2);
//   //printf("prefix2c: %.10s\n", prefix2c);
//   //printf("bprefix2: %.10s\n", bprefix2);
//   //printf("bprefix2c: %.10s\n", bprefix2c);
  
  // set up the containers for the offsets
  int p1trueForwardOffsets[TEMP_OFFSETS_BUFF][2], numP1tfo;
  int p1compForwardOffsets[TEMP_OFFSETS_BUFF][2], numP1cfo;
  int p1trueBackwardOffsets[TEMP_OFFSETS_BUFF][2], numP1tbo;
  int p1compBackwardOffsets[TEMP_OFFSETS_BUFF][2], numP1cbo;
  int p2trueForwardOffsets[TEMP_OFFSETS_BUFF][2], numP2tfo;
  int p2compForwardOffsets[TEMP_OFFSETS_BUFF][2], numP2cfo;
  int p2trueBackwardOffsets[TEMP_OFFSETS_BUFF][2], numP2tbo;
  int p2compBackwardOffsets[TEMP_OFFSETS_BUFF][2], numP2cbo;
  
  // find the offsets
  numP1tfo = OutputIndicies(pRoot, prefix1, kmer, p1trueForwardOffsets);
  numP1cfo = OutputIndicies(pRoot, prefix1c, kmer, p1compForwardOffsets);
  numP1tbo = OutputIndicies(pRoot, bprefix1, kmer, p1trueBackwardOffsets);
  numP1cbo = OutputIndicies(pRoot, bprefix1c, kmer, p1compBackwardOffsets);
  numP2tfo = OutputIndicies(pRoot, prefix2, kmer, p2trueForwardOffsets);
  numP2cfo = OutputIndicies(pRoot, prefix2c, kmer, p2compForwardOffsets);
  numP2tbo = OutputIndicies(pRoot, bprefix2, kmer, p2trueBackwardOffsets);
  numP2cbo = OutputIndicies(pRoot, bprefix2c, kmer, p2compBackwardOffsets);
  
  
  // find the potential placements
  placement_t tempPlacement;
  int insertLength = 0;
  
  double filter1, filter2, normalLikelihood;
  for(i = 0; i < numP1tfo; i++){
    filter1 = PlacementLikelihoodFTQ(theRead.readSubSeq[0], theAssembly.assemblyParts[p1trueForwardOffsets[i][0]], p1trueForwardOffsets[i][1]); // find the likelihood of the placement
    //printf("Comparing %i:%i -> %f.\n", i, p1trueForwardOffsets[i][1], filter1);
    if(filter1 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
      
      if(libProbs.IOs[0].prob > 0){ // look for BTQ placement within some standard deviation of the filter area -> <-
    for(j = 0; j < numP2tbo; j++){ // look at all of the offsets
      //printf("Comparing %i:%i -> %f and %i:%i.\n", i, p1trueForwardOffsets[i][1], filter1, j , p2trueBackwardOffsets[j][1]);
      if(p1trueForwardOffsets[i][0] == p2trueBackwardOffsets[j][0]){ // make sure they are on the same assembly part
        insertLength = p2trueBackwardOffsets[j][1] - p1trueForwardOffsets[i][1] + kmer - 1;
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[0].sigma, (float)insertLength - libProbs.IOs[0].mu); // calculate insert likelihood
        if(normalLikelihood > LIKELIHOOD_THRESHOLD){ // make sure they are within the cutoff area
          //printf("Passed normal test with %f.\n", normalLikelihood);
          filter2 = PlacementLikelihoodBTQ(theRead.readSubSeq[1], theAssembly.assemblyParts[p2trueBackwardOffsets[j][0]], p2trueBackwardOffsets[j][1] + kmer - 1); // find the likelihood of the placement
          //printf("Comparing %i:%i -> %f and %i:%i -> %f.\n", i, p1trueForwardOffsets[i][1], filter1, j , p2trueBackwardOffsets[j][1], filter2);
          if(filter2 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
        
        if(filter1*filter2*libProbs.IOs[0].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = p1trueForwardOffsets[i][1];
          tempPlacement.offset2 = p2trueBackwardOffsets[j][1];
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[0].prob*normalLikelihood;
          tempPlacement.placeInfo = 1 + 4 + 8; // oriented -> <- and true sequence and read 1 first
          tempPlacement.assemPart = p1trueForwardOffsets[i][0];
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        }
          }
        }
      }
    }
      }
      if(libProbs.IOs[1].prob > 0){ // look for BTQ placement within some standard deviation of the filter area <- ->
    for(j = 0; j < numP2tbo; j++){ // look at all of the offsets
      //printf("Comparing %i:%i -> %f and %i:%i.\n", i, p1trueForwardOffsets[i][1], filter1, j , p2trueBackwardOffsets[j][1]);
      if(p1trueForwardOffsets[i][0] == p2trueBackwardOffsets[j][0]){ // make sure they are on the same assembly part
        insertLength = p1trueForwardOffsets[i][1] + theRead.readSubSeq[0].seqLen - (p2trueBackwardOffsets[j][1] - theRead.readSubSeq[1].seqLen + kmer - 1);
        //printf("insert length: %i\n", insertLength);
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[1].sigma, (float)insertLength - libProbs.IOs[1].mu); // calculate insert likelihood
        if(normalLikelihood > LIKELIHOOD_THRESHOLD){ // make sure they are within the cutoff area
          //printf("Passed normal test!4! with %f.\n", normalLikelihood);
          filter2 = PlacementLikelihoodBTQ(theRead.readSubSeq[1], theAssembly.assemblyParts[p2trueBackwardOffsets[j][0]], p2trueBackwardOffsets[j][1] + kmer - 1); // find the likelihood of the placement
          //printf("Comparing %i:%i -> %f and %i:%i -> %f.\n", i, p1trueForwardOffsets[i][1], filter1, j , p2trueBackwardOffsets[j][1], filter2);
          if(filter2 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
        
        if(filter1*filter2*libProbs.IOs[1].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = p1trueForwardOffsets[i][1];
          tempPlacement.offset2 = p2trueBackwardOffsets[j][1];
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[1].prob*normalLikelihood;
          tempPlacement.placeInfo = 0 + 4 + 0; // oriented <- -> and true sequence and read 2 first
          tempPlacement.assemPart = p1trueForwardOffsets[i][0];
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        }
          }
        }
      }
    }
      }
    } 
  }
  
  for(i = 0; i < numP1cfo; i++){
    filter1 = PlacementLikelihoodFCQ(theRead.readSubSeq[0], theAssembly.assemblyParts[p1compForwardOffsets[i][0]], p1compForwardOffsets[i][1]); // find the likelihood of the placement
    //printf("Comparing %i:%i -> %f.\n", i, p1compForwardOffsets[i][1], filter1);
    if(filter1 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
      
      if(libProbs.IOs[0].prob > 0){ // look for BTQ placement within some standard deviation of the filter area -> <-
    for(j = 0; j < numP2cbo; j++){ // look at all of the offsets
      if(p1compForwardOffsets[i][0] == p2compBackwardOffsets[j][0]){ // make sure they are on the same assembly part
        //printf("Comparing %i:%i -> %f and %i:%i.\n", i, p1compForwardOffsets[i][1], filter1, j , p2compBackwardOffsets[j][1]);
        insertLength = p2compBackwardOffsets[j][1] + kmer - 1 - p1compForwardOffsets[i][1];
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[0].sigma, (float)insertLength - libProbs.IOs[0].mu); // calculate insert likelihood
        if(normalLikelihood > LIKELIHOOD_THRESHOLD){ // make sure they are within the cutoff area
          //printf("Passed normal test with %f.\n", normalLikelihood);
          filter2 = PlacementLikelihoodBCQ(theRead.readSubSeq[1], theAssembly.assemblyParts[p2compBackwardOffsets[j][0]], p2compBackwardOffsets[j][1] + kmer - 1); // find the likelihood of the placement
          //printf("Comparing %i:%i -> %f and %i:%i -> %f.\n", i, p1compForwardOffsets[i][1], filter1, j , p2compBackwardOffsets[j][1], filter2);
          if(filter2 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
        
        if(filter1*filter2*libProbs.IOs[0].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = p1compForwardOffsets[i][1];
          tempPlacement.offset2 = p2compBackwardOffsets[j][1];
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[0].prob*normalLikelihood;
          tempPlacement.placeInfo = 1 + 0 + 8; // oriented -> <- and comp sequence and read 1 first
          tempPlacement.assemPart = p1compForwardOffsets[i][0];
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        }
          }
        }
      }
    }
      }
      if(libProbs.IOs[1].prob > 0){ // look for BTQ placement within some standard deviation of the filter area <- ->
    for(j = 0; j < numP2cbo; j++){ // look at all of the offsets
      if(p1compForwardOffsets[i][0] == p2compBackwardOffsets[j][0]){ // make sure they are on the same assembly part
        //printf("Comparing %i:%i -> %f and %i:%i.\n", i, p1compForwardOffsets[i][1], filter1, j , p2compBackwardOffsets[j][1]);
        insertLength = p1compForwardOffsets[i][1] + theRead.readSubSeq[0].seqLen - (p2compBackwardOffsets[j][1] + kmer - 1 - theRead.readSubSeq[1].seqLen);
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[1].sigma, (float)insertLength - libProbs.IOs[1].mu); // calculate insert likelihood
        if(normalLikelihood > LIKELIHOOD_THRESHOLD){ // make sure they are within the cutoff area
          //printf("Passed normal test with %f.\n", normalLikelihood);
          filter2 = PlacementLikelihoodBCQ(theRead.readSubSeq[1], theAssembly.assemblyParts[p2compBackwardOffsets[j][0]], p2compBackwardOffsets[j][1] + kmer - 1); // find the likelihood of the placement
          //printf("Comparing %i:%i -> %f and %i:%i -> %f.\n", i, p1compForwardOffsets[i][1], filter1, j , p2compBackwardOffsets[j][1], filter2);
          if(filter2 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
        
        if(filter1*filter2*libProbs.IOs[1].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = p1compForwardOffsets[i][1];
          tempPlacement.offset2 = p2compBackwardOffsets[j][1];
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[1].prob*normalLikelihood;
          tempPlacement.placeInfo = 0 + 0 + 0; // oriented <- -> and comp sequence and read 2 first
          tempPlacement.assemPart = p1compForwardOffsets[i][0];
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        }
          }
        }
      }
    }
      }
    } 
  }
  
  for(i = 0; i < numP1tbo; i++){
    filter1 = PlacementLikelihoodBTQ(theRead.readSubSeq[0], theAssembly.assemblyParts[p1trueBackwardOffsets[i][0]], p1trueBackwardOffsets[i][1] + kmer - 1); // find the likelihood of the placement
    //printf("Comparing %i:%i -> %f.\n", i, p1trueBackwardOffsets[i][1], filter1);
    if(filter1 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
      
      if(libProbs.IOs[0].prob > 0){ // look for BTQ placement within some standard deviation of the filter area -> <-
    for(j = 0; j < numP2tfo; j++){ // look at all of the offsets
      if(p1trueBackwardOffsets[i][0] == p2trueForwardOffsets[j][0]){ // make sure they are on the same assembly part
        //printf("Comparing %i:%i -> %f and %i:%i.\n", i, p1trueBackwardOffsets[i][1], filter1, j , p2trueForwardOffsets[j][1]);
        insertLength = p1trueBackwardOffsets[i][1] + kmer - 1 - p2trueForwardOffsets[j][1];
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[0].sigma, (float)insertLength - libProbs.IOs[0].mu); // calculate insert likelihood
        if(normalLikelihood > LIKELIHOOD_THRESHOLD){ // make sure they are within the cutoff area
          //printf("Passed normal test with %f.\n", normalLikelihood);
          filter2 = PlacementLikelihoodFTQ(theRead.readSubSeq[1], theAssembly.assemblyParts[p2trueForwardOffsets[j][0]], p2trueForwardOffsets[j][1]); // find the likelihood of the placement
          //printf("Comparing %i:%i -> %f and %i:%i -> %f.\n", i, p1trueBackwardOffsets[i][1], filter1, j , p2trueForwardOffsets[j][1], filter2);
          if(filter2 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
        
        if(filter1*filter2*libProbs.IOs[0].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = p1trueBackwardOffsets[i][1];
          tempPlacement.offset2 = p2trueForwardOffsets[j][1];
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[0].prob*normalLikelihood;
          tempPlacement.placeInfo = 0 + 4 + 8; // oriented -> <- and true sequence read 2 first
          tempPlacement.assemPart = p1trueBackwardOffsets[i][0];
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        }
          }
        }
      }
    }
      }
      if(libProbs.IOs[1].prob > 0){ // look for BTQ placement within some standard deviation of the filter area <- ->
    for(j = 0; j < numP2tfo; j++){ // look at all of the offsets
      if(p1trueBackwardOffsets[i][0] == p2trueForwardOffsets[j][0]){ // make sure they are on the same assembly part
        //printf("Comparing %i:%i -> %f and %i:%i.\n", i, p1trueBackwardOffsets[i][1], filter1, j , p2trueForwardOffsets[j][1]);
        insertLength = p2trueForwardOffsets[j][1] + theRead.readSubSeq[1].seqLen - (p1trueBackwardOffsets[i][1] + kmer - 1 - theRead.readSubSeq[0].seqLen);
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[1].sigma, (float)insertLength - libProbs.IOs[1].mu); // calculate insert likelihood
        if(normalLikelihood > LIKELIHOOD_THRESHOLD){ // make sure they are within the cutoff area
          //printf("Passed normal test with %f.\n", normalLikelihood);
          filter2 = PlacementLikelihoodFTQ(theRead.readSubSeq[1], theAssembly.assemblyParts[p2trueForwardOffsets[j][0]], p2trueForwardOffsets[j][1]); // find the likelihood of the placement
          //printf("Comparing %i:%i -> %f and %i:%i -> %f.\n", i, p1trueBackwardOffsets[i][1], filter1, j , p2trueForwardOffsets[j][1], filter2);
          if(filter2 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
        if(filter1*filter2*libProbs.IOs[1].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = p1trueBackwardOffsets[i][1];
          tempPlacement.offset2 = p2trueForwardOffsets[j][1];
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[1].prob*normalLikelihood;
          tempPlacement.placeInfo = 1 + 4 + 0; // oriented <- -> and true sequence and read 1 first
          tempPlacement.assemPart = p1trueBackwardOffsets[i][0];
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        }
          }
        }
      }
    }
      }
    } 
  }
  
  for(i = 0; i < numP1cbo; i++){
    filter1 = PlacementLikelihoodBCQ(theRead.readSubSeq[0], theAssembly.assemblyParts[p1compBackwardOffsets[i][0]], p1compBackwardOffsets[i][1] + kmer - 1); // find the likelihood of the placement
    //printf("Comparing %i:%i -> %f.\n", i, p1compBackwardOffsets[i][1], filter1);
    if(filter1 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
      
      if(libProbs.IOs[0].prob > 0){ // look for BTQ placement within some standard deviation of the filter area -> <-
    for(j = 0; j < numP2cfo; j++){ // look at all of the offsets
      if(p1compBackwardOffsets[i][0] == p2compForwardOffsets[j][0]){ // make sure they are on the same assembly part
        //printf("Comparing %i:%i -> %f and %i:%i.\n", i, p1compBackwardOffsets[i][1], filter1, j , p2compForwardOffsets[j][1]);
        insertLength = p1compBackwardOffsets[i][1] + kmer - 1 - p2compForwardOffsets[j][1];
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[0].sigma, (float)insertLength - libProbs.IOs[0].mu); // calculate insert likelihood
        if(normalLikelihood > LIKELIHOOD_THRESHOLD){ // make sure they are within the cutoff area
          //printf("Passed normal test with %f.\n", normalLikelihood);
          filter2 = PlacementLikelihoodFCQ(theRead.readSubSeq[1], theAssembly.assemblyParts[p2compForwardOffsets[j][0]], p2compForwardOffsets[j][1]); // find the likelihood of the placement
          //printf("Comparing %i:%i -> %f and %i:%i -> %f.\n", i, p1compBackwardOffsets[i][1], filter1, j , p2compForwardOffsets[j][1], filter2);
          if(filter2 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
        
        if(filter1*filter2*libProbs.IOs[0].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = p1compBackwardOffsets[i][1];
          tempPlacement.offset2 = p2compForwardOffsets[j][1];
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[0].prob*normalLikelihood;
          tempPlacement.placeInfo = 0 + 0 + 8; // oriented -> <- and comp sequence and read 2 first
          tempPlacement.assemPart = p1compBackwardOffsets[i][0];
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        }
          }
        }
      }
    }
      }
      
      if(libProbs.IOs[1].prob > 0){ // look for BTQ placement within some standard deviation of the filter area <- ->
    for(j = 0; j < numP2cfo; j++){ // look at all of the offsets
      if(p1compBackwardOffsets[i][0] == p2compForwardOffsets[j][0]){ // make sure they are on the same assembly part
        //printf("Comparing %i:%i -> %f and %i:%i.\n", i, p1compBackwardOffsets[i][1], filter1, j , p2compForwardOffsets[j][1]);
        insertLength = p2compForwardOffsets[j][1] + theRead.readSubSeq[1].seqLen - (p1compBackwardOffsets[i][1] + kmer - 1 - theRead.readSubSeq[0].seqLen);
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[1].sigma, (float)insertLength - libProbs.IOs[1].mu); // calculate insert likelihood
        if(normalLikelihood > LIKELIHOOD_THRESHOLD){ // make sure they are within the cutoff area
          //printf("Passed normal test with %f.\n", normalLikelihood);
          filter2 = PlacementLikelihoodFCQ(theRead.readSubSeq[1], theAssembly.assemblyParts[p2compForwardOffsets[j][0]], p2compForwardOffsets[j][1]); // find the likelihood of the placement
          //printf("Comparing %i:%i -> %f and %i:%i -> %f.\n", i, p1compBackwardOffsets[i][1], filter1, j , p2compForwardOffsets[j][1], filter2);
          if(filter2 > LIKELIHOOD_THRESHOLD){ // see if it meets the threshold
        
        if(filter1*filter2*libProbs.IOs[1].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = p1compBackwardOffsets[i][1];
          tempPlacement.offset2 = p2compForwardOffsets[j][1];
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[1].prob*normalLikelihood;
          tempPlacement.placeInfo = 1 + 0 + 0; // oriented <- -> and comp sequence and read 1 first
          tempPlacement.assemPart = p1compBackwardOffsets[i][0];
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        }
          }
        }
      }
    }
      }
    } 
  }
  return theRead;
}

pairedRead_t FindBestPlacements(pairedRead_t theRead, assemblyPart_t thePart, libInsOrProbs_t libProbs){
  //return 0;
  // start by looking at the first subpart of the read at minimal overlap
  // look for possible (thresholded) placements
  int offset, offset2, i;
  double filter1, filter2, filter3;
  unsigned char filter1type;
  int thresholdPass = 0;
  placement_t tempPlacement;
  
  int offset2min, offset2max;
  float normalLikelihood;
  
  // look at all possible offsets
  
  for(offset = KMER_LENGTH - (int)theRead.readSubSeq[0].seqLen; offset < (int)thePart.seqLen - KMER_LENGTH; offset++){
    //printf("offset: %i\n", offset);
    filter1type = 0;
    thresholdPass = 0;
    filter1 = PlacementLikelihoodFTQ(theRead.readSubSeq[0], thePart, offset);
    //printf("filter1: %f\n", filter1);
    if(filter1 > LIKELIHOOD_THRESHOLD){
      filter1type = 1;
    }else{
      filter1 = PlacementLikelihoodFCQ(theRead.readSubSeq[0], thePart, offset);
      if(filter1 > LIKELIHOOD_THRESHOLD){
    filter1type = 2;
      }else{
    filter1 = PlacementLikelihoodBTQ(theRead.readSubSeq[0], thePart, offset);
    if(filter1 > LIKELIHOOD_THRESHOLD){
      filter1type = 3;
    }else{
      filter1 = PlacementLikelihoodBCQ(theRead.readSubSeq[0], thePart, offset);
      if(filter1 > LIKELIHOOD_THRESHOLD){
        filter1type = 4;
      }
    }
      }
    }
    
//     if(filter1type != 0){
//       printf("filter1type: %i\n", filter1type);
//     }
    
    if(filter1type == 1){
      if(libProbs.IOs[0].prob > 0){ // look for BTQ placement within some standard deviation of the filter area -> <-
    offset2min = intMax(offset - theRead.readSubSeq[1].seqLen + libProbs.IOs[0].mu - INS_OR_STD_THRESHOLD*libProbs.IOs[0].sigma, KMER_LENGTH - theRead.readSubSeq[1].seqLen);
    offset2max = intMin(offset - theRead.readSubSeq[1].seqLen + libProbs.IOs[0].mu + INS_OR_STD_THRESHOLD*libProbs.IOs[0].sigma, thePart.seqLen - KMER_LENGTH);
    for(offset2 = offset2min; offset2 < offset2max; offset2++){
      filter2 = PlacementLikelihoodBTQ(theRead.readSubSeq[1], thePart, offset2); // try to place the second part
      if(filter2 > 0.0){ // quick test to see if a reasonable placement was found
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[0].mu); // calculate insert likelihood
        if(filter1*filter2*libProbs.IOs[0].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = offset;
          tempPlacement.offset2 = offset2;
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[0].prob*normalLikelihood;
          tempPlacement.placeInfo = 1 + 4; // oriented -> <- and true sequence
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
          //PrintPlacements(theRead);
          //printf("1.1 Offset: %i, Offset2: %i, BTQ, %f, %f, %f, %f, %f\n",offset, offset2, filter1, filter2, libProbs.IOs[0].prob, GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[0].mu), filter1*filter2*libProbs.IOs[0].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[0].mu));
        }
      }
    }
      }
      if(libProbs.IOs[1].prob > 0){ // look for BTQ placement within some standard deviation of the filter area <- ->
    offset2min = intMax(offset + theRead.readSubSeq[0].seqLen - libProbs.IOs[1].mu - INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma + theRead.readSubSeq[1].seqLen, KMER_LENGTH - theRead.readSubSeq[1].seqLen);
    offset2max = (int)(offset + theRead.readSubSeq[0].seqLen - libProbs.IOs[1].mu + INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma + theRead.readSubSeq[1].seqLen);
    for(offset2 = offset2min; offset2 < offset2max; offset2++){
      filter2 = PlacementLikelihoodBTQ(theRead.readSubSeq[1], thePart, offset2);
      if(filter2 > 0.0){
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[1].mu);
        if(filter1*filter2*libProbs.IOs[1].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = offset;
          tempPlacement.offset2 = offset2;
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[1].prob*normalLikelihood;
          tempPlacement.placeInfo = 0 + 4; // oriented <- -> and true sequence
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
          //printf("1.2 Offset: %i, Offset2: %i, BTQ, %f, %f, %f, %f, %f\n",offset, offset2, filter1, filter2, libProbs.IOs[1].prob, GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[1].mu), filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[1].mu));
        }
      }
    }
      }
    }
    
    if(filter1type == 2){
      if(libProbs.IOs[0].prob > 0){ // look for BCQ placement within some standard deviation of the filter area -> <-
    offset2min = intMax(offset - theRead.readSubSeq[1].seqLen + libProbs.IOs[0].mu - INS_OR_STD_THRESHOLD*libProbs.IOs[0].sigma, KMER_LENGTH - theRead.readSubSeq[1].seqLen);
    offset2max = intMin(offset - theRead.readSubSeq[1].seqLen + libProbs.IOs[0].mu + INS_OR_STD_THRESHOLD*libProbs.IOs[0].sigma, thePart.seqLen - KMER_LENGTH);
    for(offset2 = offset2min; offset2 < offset2max; offset2++){
      filter2 = PlacementLikelihoodBCQ(theRead.readSubSeq[1], thePart, offset2);
      if(filter2 > 0.0){
        normalLikelihood = GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[0].mu);
        if(filter1*filter2*libProbs.IOs[0].prob*normalLikelihood > LIKELIHOOD_THRESHOLD){
          // IT PASSED ALL THE THRESHOLDS!!!!
          tempPlacement.offset1 = offset;
          tempPlacement.offset2 = offset2;
          tempPlacement.likelihood = filter1*filter2*libProbs.IOs[0].prob*normalLikelihood;
          tempPlacement.placeInfo = 1; // oriented -> <- and comp sequence
          theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
          //printf("2.1 Offset: %i, Offset2: %i, BCQ, %f, %f, %f, %f, %f\n",offset, offset2, filter1, filter2, libProbs.IOs[0].prob, GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[0].mu), filter1*filter2*libProbs.IOs[0].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[0].mu));
        }
      }
    }
      }
      if(libProbs.IOs[1].prob > 0){ // look for BCQ placement within some standard deviation of the filter area <- ->
    offset2min = intMax(offset + theRead.readSubSeq[0].seqLen - libProbs.IOs[1].mu - INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma + theRead.readSubSeq[1].seqLen, KMER_LENGTH - theRead.readSubSeq[1].seqLen);
    offset2max = (int)(offset + theRead.readSubSeq[0].seqLen - libProbs.IOs[1].mu + INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma + theRead.readSubSeq[1].seqLen);
    for(offset2 = offset2min; offset2 < offset2max; offset2++){
      filter2 = PlacementLikelihoodBCQ(theRead.readSubSeq[1], thePart, offset2);
      if(filter2 > 0.0){
      if(filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[1].mu) > LIKELIHOOD_THRESHOLD){
        // IT PASSED ALL THE THRESHOLDS!!!!
        tempPlacement.offset1 = offset;
        tempPlacement.offset2 = offset2;
        tempPlacement.likelihood = filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[1].mu);
        tempPlacement.placeInfo = 0; // oriented <- -> and comp sequence
        theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        //printf("2.2 Offset: %i, Offset2: %i, BCQ, %f, %f, %f, %f, %f\n",offset, offset2, filter1, filter2, libProbs.IOs[1].prob, GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[1].mu), filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[1].mu));
      }}
    }
      }
    }
    
    if(filter1type == 3){
      if(libProbs.IOs[0].prob > 0){ // look for FTQ placement within some standard deviation of the filter area -> <-
    offset2min = intMax(offset + theRead.readSubSeq[0].seqLen - libProbs.IOs[0].mu - INS_OR_STD_THRESHOLD*libProbs.IOs[0].sigma, KMER_LENGTH - theRead.readSubSeq[1].seqLen);
    offset2max = intMin(offset + theRead.readSubSeq[0].seqLen - libProbs.IOs[0].mu + INS_OR_STD_THRESHOLD*libProbs.IOs[0].sigma, thePart.seqLen - KMER_LENGTH);
    for(offset2 = offset2min; offset2 < offset2max; offset2++){
      filter2 = PlacementLikelihoodFTQ(theRead.readSubSeq[1], thePart, offset2);
      if(filter2 > 0.0){
      if(filter1*filter2*libProbs.IOs[0].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[0].mu) > LIKELIHOOD_THRESHOLD){
        // IT PASSED ALL THE THRESHOLDS!!!!
        tempPlacement.offset1 = offset;
        tempPlacement.offset2 = offset2;
        tempPlacement.likelihood = filter1*filter2*libProbs.IOs[0].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[0].mu);
        tempPlacement.placeInfo = 1 + 4; // oriented -> <- and true sequence
        theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        //printf("3.1 Offset: %i, Offset2: %i, FTQ, %f, %f, %f, %f, %f\n",offset, offset2, filter1, filter2, libProbs.IOs[0].prob, GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[0].mu), filter1*filter2*libProbs.IOs[0].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[0].mu));
      }}
    }
      }
      if(libProbs.IOs[1].prob > 0){ // look for FTQ placement within some standard deviation of the filter area <- ->
    //printf("SANITY CHECK, %i, %i, %i\n", offset - libProbs.IOs[1].mu - INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma, intMin(offset - libProbs.IOs[1].mu + INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma, thePart.seqLen - KMER_LENGTH), 1);
    offset2min = intMin(offset - theRead.readSubSeq[1].seqLen + libProbs.IOs[1].mu - INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma, thePart.seqLen - KMER_LENGTH);
    offset2max = intMin(offset - theRead.readSubSeq[1].seqLen + libProbs.IOs[1].mu + INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma, thePart.seqLen - KMER_LENGTH);
    for(offset2 = offset2min; offset2 < offset2max; offset2++){
      filter2 = PlacementLikelihoodFTQ(theRead.readSubSeq[1], thePart, offset2);
      if(filter2 > 0.0){
      if(filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[1].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[1].mu) > LIKELIHOOD_THRESHOLD){
        // IT PASSED ALL THE THRESHOLDS!!!!
        tempPlacement.offset1 = offset;
        tempPlacement.offset2 = offset2;
        tempPlacement.likelihood = filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[1].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[1].mu);
        tempPlacement.placeInfo = 0 + 4; // oriented <- -> and true sequence
        theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        //printf("3.2 Offset: %i, Offset2: %i, FTQ, %f, %f, %f, %f, %f\n",offset, offset2, filter1, filter2, libProbs.IOs[1].prob, GetInsertProbNormal(libProbs.IOs[1].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[1].mu), filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[1].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[1].mu));
      }}
    }
      }
    }
    
    if(filter1type == 4){
      if(libProbs.IOs[0].prob > 0){ // look for FCQ placement within some standard deviation of the filter area -> <-
    offset2min = intMax(offset + theRead.readSubSeq[0].seqLen - libProbs.IOs[0].mu - INS_OR_STD_THRESHOLD*libProbs.IOs[0].sigma, KMER_LENGTH - theRead.readSubSeq[1].seqLen);
    offset2max = (int)(offset + theRead.readSubSeq[0].seqLen - libProbs.IOs[0].mu + INS_OR_STD_THRESHOLD*libProbs.IOs[0].sigma);
    for(offset2 = offset2min; offset2 < offset2max; offset2++){
      filter2 = PlacementLikelihoodFCQ(theRead.readSubSeq[1], thePart, offset2);
      if(filter2 > 0.0){
      if(filter1*filter2*libProbs.IOs[0].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[0].mu) > LIKELIHOOD_THRESHOLD){
        // IT PASSED ALL THE THRESHOLDS!!!!
        tempPlacement.offset1 = offset;
        tempPlacement.offset2 = offset2;
        tempPlacement.likelihood = filter1*filter2*libProbs.IOs[0].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[0].mu);
        tempPlacement.placeInfo = 1 + 0; // oriented -> <- and comp sequence
        theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        //printf("4.1 Offset: %i, Offset2: %i, FCQ, %f, %f, %f, %f, %f\n",offset, offset2, filter1, filter2, libProbs.IOs[0].prob, GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[0].mu), filter1*filter2*libProbs.IOs[0].prob*GetInsertProbNormal(libProbs.IOs[0].sigma, (float)offset + (float)theRead.readSubSeq[0].seqLen - (float)offset2 - libProbs.IOs[0].mu));
      }}
    }
      }
      if(libProbs.IOs[1].prob > 0){ // look for FCQ placement within some standard deviation of the filter area <- ->
    offset2min = intMin(offset + libProbs.IOs[1].mu - INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma - theRead.readSubSeq[1].seqLen, thePart.seqLen - KMER_LENGTH);
    offset2max = (int)(offset + libProbs.IOs[1].mu + INS_OR_STD_THRESHOLD*libProbs.IOs[1].sigma - theRead.readSubSeq[1].seqLen);
    for(offset2 = offset2min; offset2 < offset2max; offset2++){
      filter2 = PlacementLikelihoodFCQ(theRead.readSubSeq[1], thePart, offset2);
      if(filter2 > 0.0){
      if(filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[1].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[1].mu) > LIKELIHOOD_THRESHOLD){
        // IT PASSED ALL THE THRESHOLDS!!!!
        tempPlacement.offset1 = offset;
        tempPlacement.offset2 = offset2;
        tempPlacement.likelihood = filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[1].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[1].mu);
        tempPlacement.placeInfo = 0 + 0; // oriented <- -> and comp sequence
        theRead.numPlacements = SuggestPlacement(theRead.placements, tempPlacement, theRead.numPlacements);
        //printf("4.2 Offset: %i, Offset2: %i, FCQ, %f, %f, %f, %f, %f\n",offset, offset2, filter1, filter2, libProbs.IOs[1].prob, GetInsertProbNormal(libProbs.IOs[1].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[1].mu), filter1*filter2*libProbs.IOs[1].prob*GetInsertProbNormal(libProbs.IOs[1].sigma, (float)offset2 + (float)theRead.readSubSeq[1].seqLen - (float)offset - libProbs.IOs[1].mu));
      }}
    }
      }
    }
  }
  return theRead;
}

int SuggestPlacement(placement_t oldPlacements[], placement_t newPlacement, unsigned char numPlacements){
  int i = 0, j = 0;
  int tempOffset1, tempOffset2;
  double tempLikelihood;
  unsigned char tempPlaceInfo;
  int tempAssemPart;
  if(numPlacements == 0){
    oldPlacements[0].offset1 = newPlacement.offset1;
    oldPlacements[0].offset2 = newPlacement.offset2;
    oldPlacements[0].likelihood = newPlacement.likelihood;
    oldPlacements[0].placeInfo = newPlacement.placeInfo;
    oldPlacements[0].assemPart = newPlacement.assemPart;
    return 1;
  }else{
    while(i < numPlacements){
      if(oldPlacements[i].likelihood < newPlacement.likelihood){ 
    tempOffset1 = oldPlacements[i].offset1;
    tempOffset2 = oldPlacements[i].offset2;
    tempLikelihood = oldPlacements[i].likelihood;
    tempPlaceInfo = oldPlacements[i].placeInfo;
    tempAssemPart = oldPlacements[i].assemPart;
    oldPlacements[i].offset1 = newPlacement.offset1;
    oldPlacements[i].offset2 = newPlacement.offset2;
    oldPlacements[i].likelihood = newPlacement.likelihood;
    oldPlacements[i].placeInfo = newPlacement.placeInfo;
    oldPlacements[i].assemPart = newPlacement.assemPart;
    if(j == 0){if(numPlacements < N_PLACEMENTS){numPlacements += 1;}}
    j = 1;
    newPlacement.offset1 = tempOffset1;
    newPlacement.offset2 = tempOffset2;
    newPlacement.likelihood = tempLikelihood;
    newPlacement.placeInfo = tempPlaceInfo;
    newPlacement.assemPart = tempAssemPart;
      }
      i++;
    }
    if(numPlacements < N_PLACEMENTS){
      oldPlacements[numPlacements].offset1 = newPlacement.offset1;
      oldPlacements[numPlacements].offset2 = newPlacement.offset2;
      oldPlacements[numPlacements].likelihood = newPlacement.likelihood;
      oldPlacements[numPlacements].placeInfo = newPlacement.placeInfo;
      oldPlacements[numPlacements].assemPart = newPlacement.assemPart;
      numPlacements++;
    }
  }
  printf("Suggested a placement, there are now %i.\n", numPlacements);
  return numPlacements;
}

double GetInsertProbNormal(const double sigma, const double point){
  //printf("Point: %f, p1: %f, p2: %f\n", point, erf((point + 0.5)/sqrt(2*sigma*sigma)), erf((point - 0.5)/sqrt(2*sigma*sigma)));
  return 0.5*(erf((abs(point) + 0.5)/sqrt(2*sigma*sigma)) - erf((abs(point) - 0.5)/sqrt(2*sigma*sigma)));
}

// try to place the read on the part with given offset forward with the true read sequence with Q values
double PlacementLikelihoodFTQ(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 1.0; // likelihood starts at 1.0
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  if(offset < 0){
    overlap = read.seqLen + offset;
  }else{
    overlap = intMin(part.seqLen - offset, read.seqLen);
  }
  if(offset < 0){
    for(index = 0; index < overlap; index++){
      // this next part could be optimized
      expErrors = expErrors*(1.0-GetExpectedError(getCharFromSeqByLoc(read.sequence, -offset + index), part.sequence[index], getQualityP(read.qval,-offset + index)));
      if(expErrors < LIKELIHOOD_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }else{
    for(index = 0; index < overlap; index++){
      expErrors = expErrors*(1.0-GetExpectedError(getCharFromSeqByLoc(read.sequence, index), part.sequence[offset + index], getQualityP(read.qval,index)));
      if(expErrors < LIKELIHOOD_THRESHOLD){
    return 0.0;
      }
    }
    return expErrors;
  }
  return 0.0;
}

// try to place the read on the part with given offset forward with the true read sequence with Q values
double PlacementLikelihoodFTQold1(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 0.0; // expected number of errors in this placement
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  if(offset < 0){
    overlap = read.seqLen + offset;
  }else{
    overlap = intMin(part.seqLen - offset, read.seqLen);
  }
  if(offset < 0){
    for(index = 0; index < overlap; index++){
      // this next part could be optimized
      expErrors += GetExpectedError(getCharFromSeqByLoc(read.sequence, -offset + index), part.sequence[index], getQualityP(read.qval,-offset + index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }else{
    for(index = 0; index < overlap; index++){
      expErrors += GetExpectedError(getCharFromSeqByLoc(read.sequence, index), part.sequence[offset + index], getQualityP(read.qval,index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }
  return 0.0;
}

// try to place the read on the part with given offset forward with the compliment read sequence with Q values
double PlacementLikelihoodFCQ(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 1.0; // likelihood
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  if(offset < 0){
    overlap = read.seqLen + offset;
  }else{
    overlap = intMin(part.seqLen - offset, read.seqLen);
  }
  if(offset < 0){
    for(index = 0; index < overlap; index++){
      // this next part could be optimized
      expErrors = expErrors*(1.0-GetExpectedError(getComplimentRes(getCharFromSeqByLoc(read.sequence, -offset + index)), part.sequence[index], getQualityP(read.qval,-offset + index)));
      if(expErrors < LIKELIHOOD_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }else{
    for(index = 0; index < overlap; index++){
      expErrors = expErrors*(1.0-GetExpectedError(getComplimentRes(getCharFromSeqByLoc(read.sequence, index)), part.sequence[offset + index], getQualityP(read.qval,index)));
      if(expErrors < LIKELIHOOD_THRESHOLD){
    return 0.0;
      }
    }
    return expErrors;
  }
  return 0.0;
}

// try to place the read on the part with given offset forward with the compliment read sequence with Q values
double PlacementLikelihoodFCQold1(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 0.0; // expected number of errors in this placement
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  if(offset < 0){
    overlap = read.seqLen + offset;
  }else{
    overlap = intMin(part.seqLen - offset, read.seqLen);
  }
  if(offset < 0){
    for(index = 0; index < overlap; index++){
      // this next part could be optimized
      expErrors += GetExpectedError(getComplimentRes(getCharFromSeqByLoc(read.sequence, -offset + index)), part.sequence[index], getQualityP(read.qval,-offset + index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }else{
    for(index = 0; index < overlap; index++){
      expErrors += GetExpectedError(getComplimentRes(getCharFromSeqByLoc(read.sequence, index)), part.sequence[offset + index], getQualityP(read.qval,index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }
  return 0.0;
}

// try to place the read on the part with given offset backward with the true read sequence with Q values
double PlacementLikelihoodBTQ(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 1.0; // likelihood
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  overlap = intMin(offset, read.seqLen);
  for(index = 0; index < overlap; index++){
    //printf("Comparing %c : %i and %c : %i.\n", getCharFromSeqByLoc(read.sequence, index), index, part.sequence[offset - index], offset - index);
      expErrors = expErrors*(1.0 - GetExpectedError(getCharFromSeqByLoc(read.sequence, index), part.sequence[offset - index], getQualityP(read.qval,index)));
      if(expErrors < LIKELIHOOD_THRESHOLD){
    return 0.0;
      }
    }
  return expErrors;
}

// try to place the read on the part with given offset backward with the true read sequence with Q values
double PlacementLikelihoodBTQold1(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 0.0; // expected number of errors in this placement
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  overlap = intMin(offset, read.seqLen);
  for(index = 0; index < overlap; index++){
    //printf("Comparing %c : %i and %c : %i.\n", getCharFromSeqByLoc(read.sequence, index), index, part.sequence[offset - index], offset - index);
      expErrors += GetExpectedError(getCharFromSeqByLoc(read.sequence, index), part.sequence[offset - index], getQualityP(read.qval,index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
  return GetLikelihoodOfPlacementByErrors(expErrors);
}

// try to place the read on the part with given offset backward with the compliment read sequence with Q values
double PlacementLikelihoodBCQ(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 1.0; // likelihood
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  overlap = intMin(offset, read.seqLen);
  for(index = 0; index < overlap; index++){
    //printf("Comparing %c : %i and %c : %i.\n", getCharFromSeqByLoc(read.sequence, index), index, part.sequence[offset - index], offset - index);
      expErrors = expErrors*(1.0-GetExpectedError(getComplimentRes(getCharFromSeqByLoc(read.sequence, index)), part.sequence[offset - index], getQualityP(read.qval,index)));
      if(expErrors < LIKELIHOOD_THRESHOLD){
    return 0.0;
      }
    }
  return expErrors;
}

// try to place the read on the part with given offset backward with the compliment read sequence with Q values
double PlacementLikelihoodBCQold1(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 0.0; // expected number of errors in this placement
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  overlap = intMin(offset, read.seqLen);
  for(index = 0; index < overlap; index++){
    //printf("Comparing %c : %i and %c : %i.\n", getCharFromSeqByLoc(read.sequence, index), index, part.sequence[offset - index], offset - index);
      expErrors += GetExpectedError(getComplimentRes(getCharFromSeqByLoc(read.sequence, index)), part.sequence[offset - index], getQualityP(read.qval,index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
  return GetLikelihoodOfPlacementByErrors(expErrors);
}

// try to place the read on the part with given offset backward with the true read sequence with Q values
double PlacementLikelihoodBTQold(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 0.0; // expected number of errors in this placement
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  if(offset < 0){
    overlap = read.seqLen + offset;
  }else{
    overlap = intMin(part.seqLen - offset, read.seqLen);
  }
  if(offset < 0){
    for(index = 0; index < overlap; index++){
      // this next part could be optimized
      expErrors += GetExpectedError(getCharFromSeqByLoc(read.sequence, read.seqLen - 1 - (-offset + index)), part.sequence[index], getQualityP(read.qval,-offset + index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }else{
    for(index = 0; index < overlap; index++){
      expErrors += GetExpectedError(getCharFromSeqByLoc(read.sequence, read.seqLen - 1 - (index)), part.sequence[offset + index], getQualityP(read.qval,index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }
  return 0.0;
}

// try to place the read on the part with given offset backward with the compliment read sequence with Q values
double PlacementLikelihoodBCQold(const readSequence_t read, const assemblyPart_t part, const int offset){
  int index;
  double expErrors = 0.0; // expected number of errors in this placement
  if(read.seqLen > part.seqLen){
    return 0; // the read cannot be longer than the assembly piece
  }
  int overlap;
  if(offset < 0){
    overlap = read.seqLen + offset;
  }else{
    overlap = intMin(part.seqLen - offset, read.seqLen);
  }
  if(offset < 0){
    for(index = 0; index < overlap; index++){
      // this next part could be optimized
      expErrors += GetExpectedError(getComplimentRes(getCharFromSeqByLoc(read.sequence, read.seqLen - 1 - (-offset + index))), part.sequence[index], getQualityP(read.qval,-offset + index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }else{
    for(index = 0; index < overlap; index++){
      expErrors += GetExpectedError(getComplimentRes(getCharFromSeqByLoc(read.sequence, read.seqLen - 1 - (index))), part.sequence[offset + index], getQualityP(read.qval,index));
      if(expErrors > PLACEMENT_THRESHOLD){
    return 0.0;
      }
    }
    return GetLikelihoodOfPlacementByErrors(expErrors);
  }
  return 0.0;
}

double GetLikelihoodOfPlacementByErrors(const double expErrors){
  //printf("Expected errors: %f\n", expErrors);
  return exp(-1.0*expErrors);
}

double GetExpectedError(const char readRes, const char partRes, const double qval){
  return (readRes == partRes)*(1.0-qval) + (partRes == 'N')*0.75 + ((readRes != partRes)*(partRes != 'N'))*(0.66666666667 + 0.33333333333*qval);
}

double GetExpectedErrorOld(const char readRes, const char partRes, const double qval){
  if(readRes == partRes){
    return 1.0 - qval; // there is a 1-qval chance it is wrong
  }else if(partRes == 'N'){
    return 0.75; // it is unknown so there is a 3/4 chance it is wrong
  }else{
    return 0.66666666667 + 0.33333333333*qval; // see writeup, this assumes uniform across unreported bases (to be modified later to account for a GC/AT bias)
  }
}

void combinePlacements(pairedRead_t theRead, placement_t tempPlacements[]){
  // search through the placements in theRead and insert in tempPlacements when they are better
}
