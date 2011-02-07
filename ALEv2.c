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

int main(int argc, char **argv){
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
}