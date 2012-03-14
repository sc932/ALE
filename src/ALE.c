/*
 * Copyright (C) 2010,2011,2012 Scott Clark. All rights reserved.
 *
 * Developed by:
 * Scott Clark
 * Cornell University Center for Applied Mathematics
 * http://cam.cornell.edu
 * AND
 * Rob Egan
 * Department of Energy Joint Genome Institute
 * http://jgi.doe.gov
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a 
 * copy of this software and associated documentation files (the "Software"), 
 * to deal with the Software without restriction, including without limitation 
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 * and/or sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 *   1. Redistributions of source code must retain the above copyright notice, 
 *      this list of conditions and the following disclaimers.
 *   2. Redistributions in binary form must reproduce the above copyright 
 *      notice, this list of conditions and the following disclaimers in the 
 *      documentation and/or other materials provided with the distribution.
 *   3. Neither the names of Cornell University, The Joint Genome Institute, 
 *      nor the names of its contributors may be used to endorse or promote 
 *      products derived from this Software without specific prior written 
 *      permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS WITH THE SOFTWARE.
 */

// For more information on the license please see 
// The University of Illinois/NCSA Open Source License
// http://www.opensource.org/licenses/UoI-NCSA.php

// ALE.c

#include "ALE.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ALEhelpers.h"
#include "ALElike.h"
#include "geneTree.h"

int main(int argc, char **argv){
    // Handles input variations
    if (argc < 2) {
        // no input, so print usage and options
        printf(USAGE, argv[0]);
        printf("%s", SHORT_OPTIONS);
        return 0;
    }
    if (argc < 4) {
        // the input was shorter than anticipated, print long options if -h specified
        if(argc <= 1 || (argv[1][0] == '-' && argv[1][1] == 'h')){
            printf("%s", WELCOME_MSG);
            printf(USAGE, argv[0]);
            printf("\n%s", LONG_OPTIONS);
            return 0;
        }else{
            // did not get enough command line input AND no help requested
            printf(USAGE, argv[0]);
            printf("%s", SHORT_OPTIONS);
            return 0;
        }
    }
    
    int options;
    int kmerLen = 4;
    char placementOut[256];
    samfile_t *placementBam = NULL;
    *placementOut = '\0';
    libraryParametersT *libParams = NULL;
    double outlierFraction = 0.02;
    int qOff = -1;
    int printAllAleOutput = 1;
    
    if(argc > 4) { // look for command line options
        for(options = 1; options < argc - 3; options++){ // search over all options
           if(strcmp(argv[options], "-kmer") == 0){
                kmerLen = atoi(argv[options+1]);
                if(kmerLen > 20){
                    printf("-kmer option of %i not in range [2,20], set to default [4].\n", kmerLen);
                    kmerLen = 4;
                }
                options++;
            }else if(strcmp(argv[options], "-qOff") == 0){
                qOff = atoi(argv[options+1]);
                if(qOff != 64 && qOff != 33 && qOff != 0){
                    printf("-qOff option of %i not in set [33,64], will be set to 33.\n", qOff);
                    qOff = 33;  // SAM/BAM specification is for ascii - 33.
                }
                options++;
            }else if(strcmp(argv[options], "-pl") == 0){
                strcpy(placementOut, argv[options+1]);
                options++;
            }else if(strcmp(argv[options], "-pm") == 0){
                libParams = malloc(sizeof(libraryParametersT));
                importLibraryParameters(libParams, argv[options+1]);
                options++;
            }else if(strcmp(argv[options], "-nout") == 0){
                printAllAleOutput = 0;
                printf("Turned off per-base output\n");
            }else if(strcmp(argv[options], "-minLL") == 0) {
            	setMinLogLike(atof(argv[options+1]));
            	options++;
            	printf("Set minLogLike to: %0.1f\n", getMinLogLike());
            } else{
                printf("Could not find option %s\n", argv[options]);
            }
        }
    }
    
    // input and output files
    printf("BAM file: %s\n", argv[argc - 3]);
    printf("Assembly fasta file: %s\n", argv[argc - 2]);
    printf("ALE Output file: %s\n", argv[argc - 1]);

    // attempt to open the bam input file
    samfile_t *ins = openSamOrBam(argv[argc - 3]);
    
    printf("Reading in assembly...\n");
    assemblyT *theAssembly = loadAssembly(argv[argc - 2]);
    if (!validateAssemblyIsSameAsAlignment(ins->header, theAssembly)) {
        printf("Error! Assembly fasta %s does not match alignment file %s!\n", argv[argc-2], argv[argc-3]);
        exit(1);
    }
    
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

    if(libParams != NULL){
      saveLibraryParameters(libParams, argv[argc - 1]); // the ALE output file name
    }

    // calculate the insert mean/std if not given
    if(libParams == NULL){
        printf("Insert length and std not given, will be calculated from input map.\n");

        libParams = computeLibraryParameters(ins, outlierFraction, qOff);
        saveLibraryParameters(libParams, argv[argc - 1]); // the ALE output file name

        // close and re-open bam file
        samclose(ins);
        ins = openSamOrBam(argv[argc - 3]);
    }

    // set assembly avgReadSize to that of library (rounded to nearest int)
    theAssembly->avgReadSize = (double)libParams->avgReadSize;

    printf("Calculating GC content of reference over average read size\n");
    calculateGCcont(theAssembly, libParams->avgReadSize);

    // place reads and compute statistics on the assembly
    printf("Computing read placements and depths\n");
    computeReadPlacements(ins, theAssembly, libParams, placementBam);
    
    // compute statistics on assembly
    printf("Computing k-mer statistics...\n");
    computeKmerStats(theAssembly, kmerLen);
    printf("Done computing k-mer statistics.\n");
    
    printf("Computing depth statistics...\n");
    computeDepthStats(theAssembly);
    printf("Done computing depth statistics.\n");

    printf("Computing statistics on expected missing...\n");
    applyExpectedMissingLength(theAssembly);
    printf("Done computing statistics on expected missing.\n");
    
    FILE *out = fopen(argv[argc - 1], "w");
    if(out == NULL){
        printf("Error! Could not open output file: %s\n", argv[argc - 1]);
    }
    writeToOutput(theAssembly, printAllAleOutput, out);
    fclose(out);
    
    printf("Output is in file: %s\n", argv[argc - 1]);

    if (placementBam != NULL) {
        printf("Closing placement file\n");
        samclose(placementBam);
    }

    printf("Closing input file\n");
    samclose(ins);
    
    if (libParams != NULL)
        free(libParams);

    freeAssembly(theAssembly);
    return 0;
}
