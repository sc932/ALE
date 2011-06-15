// (C) 2010 Scott Clark, JGI, LBNL

// Compiling instructions
// $ cc -g -O2 ALEv2.c -o ALE -lz -lm
// $ ./ALE <parameter file>
// requires: ALE.h, kseq.h, zlib.h, geneTree.h

// bowtie -t -I 0 -X 400 --fr -a -l 10 -v 3 -e 300 -S --threads 6 --sam-nohead Ecoli_complete -1 part1_EcoliReads2800k.fna -2 part2_EcoliReads2800k.fna Ecoli_complete.map.sam

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
        printf(USAGE, argv[0]);
        printf("%s", SHORT_OPTIONS);
        return 0;
    }
    if (argc < 4) {
        if(argc <= 1 || (argv[1][0] == '-' && argv[1][1] == 'h')){
            printf("%s", WELCOME_MSG);
            printf(USAGE, argv[0]);
            printf("\n%s", LONG_OPTIONS);
            return 0;
        }else{
            printf(USAGE, argv[0]);
            printf("%s", SHORT_OPTIONS);
            return 0;
        }
    }
    
    int options;
    int numberAssemblyPieces = 0;
    int kmerLen = 4;
    double insertLength = -1.0;
    double insertStd = -1.0;
    int avgReadSize = 0;
    int qOff = 33;
    char placementOut[256];
    samfile_t *placementBam = NULL;
    *placementOut = '\0';
    double chimerFraction = 0.005;
    double outlierFraction = 0.02;
    
    if(argc > 5) { // look for command line options
        for(options = 1; options < argc - 4; options++){ // search over all options
            if(strcmp(argv[options], "-nap") == 0){
                numberAssemblyPieces = atoi(argv[options+1]);
                options++;
            }else if(strcmp(argv[options], "-kmer") == 0){
                kmerLen = atoi(argv[options+1]);
                if(kmerLen > 20){
                    printf("-kmer option of %i not in range [2,20], set to default [4].\n", kmerLen);
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
	        }else if(strcmp(argv[options], "-ars") == 0){
                avgReadSize = atoi(argv[options+1]);
                if(avgReadSize <= 0){
                    printf("-ars option of %i not in range (1,inf), will be calculated from input.\n", avgReadSize);
                    avgReadSize = -1.0;
                }
                options++;
	        }else if(strcmp(argv[options], "-qOff") == 0){
                qOff = atoi(argv[options+1]);
                if(qOff != 64 && qOff != 33){
                    printf("-qOff option of %i not in set [33,64], will be set to 33.\n", qOff);
                    qOff = 33;  // SAM/BAM specification is for ascii - 33.
                }
                options++;
	        }else if(strcmp(argv[options], "-pl") == 0){
	        	strcpy(placementOut, argv[options+1]);
                options++;
	        } else if(strcmp(argv[options], "-chi") == 0){
	        	double inChimerFraction = atof(argv[options+1]);
	        	if (inChimerFraction >= 1.0 || inChimerFraction < 0.0) {
	        		printf("-chi option of %f not in range[0,1).  Resetting to %f%%\n", inChimerFraction, chimerFraction);
	        	} else
	        		chimerFraction = inChimerFraction;
	        	options++;
            }else{
                printf("Could not find option %s\n", argv[options]);
            }
        }
    }
    
    // input and output files
    printf("BAM file: %s\n", argv[argc - 3]);
    printf("Assembly fasta file: %s\n", argv[argc - 2]);
    printf("ALE Output file: %s\n", argv[argc - 1]);

    // attempt to open the bam input file
    samfile_t *ins = samopen(argv[argc - 3], "rb", 0);
    if (ins == 0) {
    	printf("Error! Failed to open BAM file %s\n", argv[argc - 3]);
    	exit(1);
    }

    // attempt to open the input file
    gzFile *assemblyFile = gzopen(argv[argc - 2], "r");
    kseq_t *Aseq;
    if(assemblyFile == NULL){
        printf("Error! Could not open assembly file: %s\n", argv[argc - 2]);
        exit(1);
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
    
    // initialize
    int i;
    double likelihoodRead1, likelihoodRead2, likelihoodInsert;
    alignSet_t alignments[N_PLACEMENTS];
    bam1_t *samReadPairs[N_PLACEMENTS*2];
    for(i=0; i < N_PLACEMENTS; i++) {
    	samReadPairs[i*2] = bam_init1();
    	samReadPairs[i*2+1] = bam_init1();
    }
    int samReadPairIdx = 0;
    
    alignSet_t *thisAlignment = NULL;
    alignSet_t *currentAlignment = NULL;
    alignSet_t *head = currentAlignment;
    
    //printAssembly(theAssembly);
    assemblySanityCheck(theAssembly);
    
    if (*placementOut != '\0') {
        printf("Placement file: %s\n", placementOut);
        char *mode = "wbu";
        if (placementOut[0] != '-' )
        	mode[2] = '\0'; // compress bam if not part of a pipe
        placementBam = samopen(placementOut, mode, ins->header);
        if(placementBam == 0){
            printf("Error! Could not write to the placement BAM file: %s\n", placementOut);
            exit(1);
        }
    }

    printf("Reading in the map and computing statistics...\n");
    // read in the first part of the read
    int keepGoing = 1;

    // calculate the insert mean/std if not given
    
    double readSizeTotal = 0.0;
    double lengthTotal = 0.0;
    double lengthStd = 0.0;
    int readCount = 0;

    if(insertLength == -1 || insertStd == -1 || avgReadSize == 0){
        int mapLens[mapLens_MAX];
        for(i = 0; i < mapLens_MAX; i++){
            mapLens[i] = 0;
        }
//         int GCmaps[GCmaps_MAX];
//         for(i = 0; i < GCmaps_MAX; i++){
//             GCmaps[i] = 0;
//         }
        printf("Insert length or std or avg read size not given, will be calculated from input map.\n");
        while(keepGoing > 0){
        	bam1_t *thisRead = samReadPairs[samReadPairIdx*2];
        	bam1_t *thisReadMate = samReadPairs[samReadPairIdx*2+1];
        	int numReads = readMatesBAM(ins, thisRead, thisReadMate);
        	if (numReads == 0)
        		break;
        	else if (numReads == 1)
        		continue; // not a pair

            if ((thisRead->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) != 0) {
            	//printf("WARNING: %s and %s has a mate that is unmapped.\n", bam1_qname(thisRead), bam1_qname(thisReadMate));
            	continue; // at least one mate does not map
            }

            if (strcmp(bam1_qname(thisRead), bam1_qname(thisReadMate)) != 0) {
            	printf("WARNING: read invalid mate pair: %s %s\n", bam1_qname(thisRead), bam1_qname(thisReadMate));
            	continue;
            }

            readSizeTotal += getSeqLenBAM(thisRead) + getSeqLenBAM(thisReadMate);
            //GCtot = getGCtotal(read.readSeq, getSeqLen(read.readSeq), readMate.readSeq, getSeqLen(readMate.readSeq));
            int mapLen = getMapLenBAM(thisRead, thisReadMate);
            lengthTotal += mapLen;
            if (mapLen < mapLens_MAX)
                mapLens[mapLen] += 1;
            //GCmaps[GCtot] += 1;
            if ((++readCount & 0xffff) == 0)
            	printf("Read %d reads\n", readCount);
        }

        // zero out top and bottom outliers
        int observed = 0;
        int purged = 0;
        for(i = 0; i < mapLens_MAX; i++) {
        	if (observed > (1.0-outlierFraction) * readCount) {
        	    purged += mapLens[i];
        	    lengthTotal -= i*mapLens[i];
        	    mapLens[i] = 0;
        	} else {
        		observed += mapLens[i];
        	}
        	if (observed < outlierFraction * readCount) {
        		purged += mapLens[i];
        		lengthTotal -= i*mapLens[i];
        		mapLens[i] = 0;
        	}
        }
        printf("Read %i properly mated read pairs, purged %i %d%% & %d%% outliers\n", readCount, purged, (int) outlierFraction*100, (int) (1.0-outlierFraction)*100);
        int modifiedReadCount = readCount - purged;

        avgReadSize = (int)(readSizeTotal/((double)readCount*2.0));
	    printf("Found sample avg read size to be %i\n", avgReadSize);
	    if(insertLength == -1 || insertStd == -1){
	        insertLength = lengthTotal/(double)modifiedReadCount;
	        printf("Found sample avg insert length to be %f from %i mapped reads\n", insertLength, modifiedReadCount);

	        for(i = 0; i < mapLens_MAX; i++){
	            if(mapLens[i] > 0){
		            lengthStd += mapLens[i]*((double)i - insertLength)*((double)i - insertLength);
		            printf("i : mapLens[i] :: %i : %i\n", i, mapLens[i]);
	            }
	        }
	        insertStd = sqrt(lengthStd/(double)(modifiedReadCount-1));
	        printf("Found sample insert length std to be %f\n", insertStd);
	    }

	    // close and re-open bam file
	    samclose(ins);
	    ins = samopen(argv[argc - 3], "rb", 0);
	    if (ins == 0) {
	    	printf("Error! Failed to open BAM file %s\n", argv[argc - 3]);
	    	exit(1);
	    }
    }

    int failedToPlace = 0;
    int placed = 0;

    keepGoing = 1;
    readCount = 0;
    while(keepGoing > 0){
    	bam1_t *thisRead = samReadPairs[samReadPairIdx*2];
    	bam1_t *thisReadMate = samReadPairs[samReadPairIdx*2+1];
    	thisAlignment = &alignments[samReadPairIdx];
    	samReadPairIdx++;

    	likelihoodRead1 = 1.0;
    	likelihoodRead2 = 1.0;
    	likelihoodInsert = 1.0;
    	thisAlignment->likelihood = 1.0;

        int numReads = readMatesBAM(ins, thisRead, thisReadMate);
        if ((++readCount & 0xffff) == 0)
        	printf("Read %d reads\n", readCount);

    	if (numReads == 2) {
    		// two reads
            thisAlignment->start1 = -1;
            thisAlignment->end1   = -1;
            thisAlignment->start2 = -1;
            thisAlignment->end2   = -1;
            if ((thisRead->core.flag & BAM_FUNMAP) == 0) {
            	thisAlignment->likelihood *= likelihoodRead1  = getMatchLikelihoodBAM(thisRead, qOff);
                thisAlignment->start1 = thisRead->core.pos;
                thisAlignment->end1   = bam_calend(&thisRead->core, bam1_cigar(thisRead));
            } else {
            	printf("read1 unmapped %s\n", bam1_qname(thisRead));
            	likelihoodRead1 = 0.0;
            }
            if ((thisReadMate->core.flag & BAM_FUNMAP) == 0) {
            	thisAlignment->likelihood *= likelihoodRead2  = getMatchLikelihoodBAM(thisReadMate, qOff);
                thisAlignment->start2 = thisReadMate->core.pos;
                thisAlignment->end2   =  bam_calend(&thisReadMate->core, bam1_cigar(thisReadMate));
            } else {
            	printf("read2 unmapped %s\n", bam1_qname(thisReadMate));
            	likelihoodRead2 = 0.0;
            }
            // Check mapping
    		if ((thisRead->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) == 0) {
    			// both map
                if (thisRead->core.tid == thisReadMate->core.mtid) {
                	// both map to same target
                	likelihoodInsert = getInsertLikelihoodBAM(thisRead, thisReadMate, insertLength, insertStd);
                } else {
                    printf("WARNING: chimeric read mate pair %s.\n", bam1_qname(thisRead));
                    likelihoodInsert = chimerFraction; // expected chimeric rate
                	// apply placement of read2 to target2...
                	// HACK!!!
                	// TODO fix for multiple possible placements
                	alignSet_t _read2Only;
                	alignSet_t *read2Only = &_read2Only;
                	read2Only->start1 = thisAlignment->start2;
                	read2Only->end2 = thisAlignment->end2;
                	read2Only->start2 = read2Only->end2 = -1;
                	read2Only->likelihood = likelihoodRead2 * likelihoodInsert;
                	if (read2Only->likelihood > 0.0) {
                	    strcpy(read2Only->name, bam1_qname(thisReadMate));
                	    strcpy(read2Only->mapName, getTargetName(ins, thisReadMate));
                	    read2Only->nextAlignment = NULL;
                	    int winner = applyPlacement(read2Only, theAssembly);
                	    if (winner < 0) {
                		    printf("WARNING: no placement found for read2 of chimer %s!\n", read2Only->name);
                		    continue;
                	    }
                	}
                	// then treat like single read1 to target1
                    thisAlignment->likelihood = likelihoodRead1;
                    likelihoodRead2 = 0.0;  // there is no read2
                    thisAlignment->start2 = thisAlignment->end2 = -1;
                }
                thisAlignment->likelihood *= likelihoodInsert;
    		}
    		if (!(likelihoodRead1 > 0.0 || likelihoodRead2 > 0.0)) {
    			printf("WARNING: Detected unmapped read mate pair: %s %s %f %f %f\n", bam1_qname(thisRead), bam1_qname(thisReadMate), likelihoodRead1, likelihoodRead2, likelihoodInsert);
    			continue; // no alignment to evaluate
    		}

    	} else if (numReads == 1) {
    	    printf("WARNING: detected single read %s\n", bam1_qname(thisRead));
    	    thisAlignment->likelihood *= likelihoodRead1 = getMatchLikelihoodBAM(thisRead, qOff);
            thisAlignment->start1 = thisRead->core.pos;
            thisAlignment->end1   = bam_calend(&thisRead->core, bam1_cigar(thisRead));
            likelihoodRead2 = 0.0;  // there is no read2
            likelihoodInsert = 0.0; // there is no insert
            thisAlignment->start2 = thisAlignment->end2 = -1;
    	} else {
    	    break; // no more reads
    	}
    	strcpy(thisAlignment->name, bam1_qname(thisRead));
    	strcpy(thisAlignment->mapName, getTargetName(ins, thisRead));
    	thisAlignment->nextAlignment = NULL;

    	//printf("Likelihoods (%s): %12f %12f %12f\n", bam1_qname(thisRead), likelihoodRead1, likelihoodRead2, likelihoodInsert);
        //printf("%s : %s .\n", currentAlignment->name, read.readName);

        if(currentAlignment == NULL || head == NULL){ // first alignment
            //printf("First alignment.\n");
        	currentAlignment = head = thisAlignment;
        }else if(strcmp(head->name, thisAlignment->name) == 0){ // test to see if this is another alignment of the current set or a new one
            // extend the set of alignments
        	//printf("Same alignment!\n");
        	currentAlignment->nextAlignment = thisAlignment;
            currentAlignment = thisAlignment;
        }else{ // new alignment
            //printf("New alignment!\n");
            // do the statistics on *head, that read is exhausted
            // printAlignments(head);
         	int winner;
            if((winner = applyPlacement(head, theAssembly)) == -1){
                failedToPlace++;
		    } else {
		      	if (placementBam != NULL) {
		      		bam_write1(placementBam->x.bam, samReadPairs[winner*2]);
		      		bam_write1(placementBam->x.bam, samReadPairs[winner*2+1]);
		      	}
		       	placed++;
		    }
		    //printf("%s : %f\n", readMate.readName, head->likelihood);
            //printf("Winner is %d of %d for %s at %f. Next is %s\n", winner, samReadPairIdx-1, head->name, alignments[winner].likelihood, thisAlignment->name);

            // refresh head and current alignment
            currentAlignment = &alignments[0];
            copyAlignment(currentAlignment, thisAlignment);
            samReadPairIdx = 1;
            head = currentAlignment;
        }

        // make sure we do not overflow the number of placements
    	if (samReadPairIdx >= N_PLACEMENTS) {
    		//printf("WARNING: Exceeded maximum number of duplicate placements: %s\n", thisAlignment->name);
    		int previous = N_PLACEMENTS-2;
    		currentAlignment = &alignments[previous];
    		alignSet_t *tmp = head;
    		alignSet_t *leastLikely = head;
    		double least = head->likelihood;
    		while (tmp != NULL) {
    			if (tmp->likelihood < least) {
    				least = tmp->likelihood;
    				leastLikely = tmp;
    			}
    			tmp = tmp->nextAlignment;
    		}
    		if (thisAlignment->likelihood > least) {
    			// overwrite previous with current
    			printf("WARNING: exceeded maximum placements. Replacing %f with %f\n", leastLikely->likelihood, thisAlignment->likelihood);
    			tmp = leastLikely->nextAlignment;
    			copyAlignment(leastLikely, thisAlignment);
    			leastLikely->nextAlignment = tmp;
    		} else {
    			printf("WARNING: exceeded maximum placements.  Dropping low probability placement %f\n", thisAlignment->likelihood);
    		}
    		currentAlignment->nextAlignment = NULL;
    		samReadPairIdx = N_PLACEMENTS-1;
    		continue;
    	}

    }
    
    printf("%i maps placed, %i maps failed to place.\n", placed, failedToPlace);
    
//     //printAssembly(theAssembly);
//     assemblySanityCheck(theAssembly);
//     
//     //compute GC content
// 
//     
//     // clean up the final alignment
//     applyPlacement(head, theAssembly);
//     //printAlignments(head);
//     
//     calculateGCcont(theAssembly, avgReadSize);
//     FILE *outGC = fopen("aleGCcont", "w");
//     int i, j;
//     if(theAssembly->numContigs > 1){
//         for(i = 0; i < theAssembly->numContigs; i++){
//             fprintf(outGC, ">%s : length=%i\n", theAssembly->contigs[i].name, theAssembly->contigs[i].seqLen);
//             for(j = 0; j < theAssembly->contigs[i].seqLen; j++){
//                 fprintf(outGC, "%f,%f\n", theAssembly->contigs[i].depth[j], theAssembly->contigs[i].GCcont[j]);
//             }
//         }
//     }else{
//         fprintf(outGC, ">%s %i > depth ln(depthLike) ln(placeLike) ln(kmerLike) ln(totalLike)\n", theAssembly->contigs->name, theAssembly->contigs->seqLen);
//         for(j = 0; j < theAssembly->contigs->seqLen; j++){
//             fprintf(outGC, "%f,%f\n", theAssembly->contigs->depth[j], theAssembly->contigs->GCcont[j]);
//         }
//     }
//     fclose(outGC);
//     
//     printf("Done reading in the map.\n");
    
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
    fclose(out);
    
    printf("Done computing statistics.\nOutput is in file: %s\n", argv[argc - 1]);

    // tear down SAM/BAM variables
    for(i=0; i < N_PLACEMENTS; i++) {
    	if (samReadPairs[i*2] != NULL)
    	    bam_destroy1(samReadPairs[i*2]);
    	if (samReadPairs[i*2+1] != NULL)
    		bam_destroy1(samReadPairs[i*2+1]);
    }
    if (placementBam != NULL) {
        printf("Closing placement file\n");
        samclose(placementBam);
    }

    printf("Closing input file\n");
    samclose(ins);
    
    free(theAssembly);
    return 0;
}

