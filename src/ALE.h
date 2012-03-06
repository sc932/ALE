// ALE.h

#ifndef _ALE_H_
#define _ALE_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ALEhelpers.h"
#include "ALElike.h"
#include "geneTree.h"

static const char WELCOME_MSG[80] = "Welcome to the Assembly Likelihood Estimator!\n(C) 2010 Scott Clark\n\n";
static const char USAGE[80] = "Usage: %s [-options] readSorted.[s|b]am assembly.fasta[.gz] ALEoutput.txt\n";
static const char SHORT_OPTIONS[80] = "    Options:\n    -h : print out help\n";
static const char LONG_OPTIONS[1024] = "Options: <i>nt <f>loat <s>tring [default]\n  -h        : print out this help\n\n  -kmer <f> : Kmer depth for kmer stats [4]\n  -qOff <i> : Quality ascii offset (illumina) [33] or 64 (or 0)\n  -pl <s>   : placementOutputBAM\n  -pm <s>   : library parameter file (auto outputs .param)\n  -nout     : only output meta information (no per base) [off]\n  -minLL    : the minimum log Likelihood (-120)\n\n\n";

#endif
