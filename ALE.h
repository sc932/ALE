// ALE.h

#ifndef _ALE_H_
#define _ALE_H_

#include <zlib.h>
#include <assert.h>
#include "kseq.h"
#include "sam.h"
#include "bam.h"

KSEQ_INIT(gzFile, gzread);

/**********************************************
**** CONSTANTS - MUST BE SET EVERY TIME!!! ****
**********************************************/

// the lengths must be a multiple of 4 (1048576 = 2^20)
#define MAX_READ_LENGTH             256
#define MAX_LINE_LENGTH             2048
#define MAX_CONTIG_LENGTH           12800000
#define MAX_FORWARD_SUBSEEDS        10000

//#define NUM_PAIRED_READS_ON_NODE    4000
//#define NUM_ASSEMBLY_PARTS          1
#define NUM_NODES                   1
#define KMER_LENGTH                 50
#define MIN_DEPTH                   1
#define MAX_NUM_READS_APPENDS_INTO  10000
#define N_PLACEMENTS                100
#define PLACEMENT_THRESHOLD         4
#define LIKELIHOOD_THRESHOLD        0.000000001
#define INS_OR_STD_THRESHOLD        10.0
#define TEMP_OFFSETS_BUFF           400
#define POISSON_COMB                1

static unsigned char powup[8] = {1,2,4,8,16,32,64,128};

// http://en.wikipedia.org/wiki/FASTQ_format
// generated from scripts.py
static double QtoP[63] = {0.0,0.205671765276,0.36904265552,0.498812766373,0.601892829447,0.683772233983,0.748811356849,0.800473768503,0.841510680754,0.874107458821,0.9,0.920567176528,0.936904265552,0.949881276637,0.960189282945,0.968377223398,0.974881135685,0.98004737685,0.984151068075,0.987410745882,0.99,0.992056717653,0.993690426555,0.994988127664,0.996018928294,0.99683772234,0.997488113568,0.998004737685,0.998415106808,0.998741074588,0.999,0.999205671765,0.999369042656,0.999498812766,0.999601892829,0.999683772234,0.999748811357,0.999800473769,0.999841510681,0.999874107459,0.9999,0.999920567177,0.999936904266,0.999949881277,0.999960189283,0.999968377223,0.999974881136,0.999980047377,0.999984151068,0.999987410746,0.99999,0.999992056718,0.999993690427,0.999994988128,0.999996018928,0.999996837722,0.999997488114,0.999998004738,0.999998415107,0.999998741075,0.999999,0.999999205672,0.999999369043};

const static double lnfactconst = 0.918938533204672741780329;

static const char WELCOME_MSG[80] = "Welcome to the Assembly Likelihood Estimator!\n(C) 2010 Scott Clark\n\n";
static const char USAGE[80] = "Usage: %s [-options] readSorted.[s|b]am assembly.fasta[.gz] ALEoutput.txt\n";
static const char SHORT_OPTIONS[80] = "    Options:\n    -h : print out help\n";
static const char LONG_OPTIONS[1024] = "Options: <i>nt <f>loat [default]\n  -h        : print out this help\n\n   -kmer <f> : Kmer depth for kmer stats [4]\n  -qOff <i> : Quality ascii offset (illumina) [33] or 64 (or 0)\n  -pl placementOutputBAM\n\n\n";

/*****************************
**** STRUCTS FOR THE TREE ****
*****************************/

struct IOtriplet{
  float prob;
  float mu;
  float sigma;
};

struct libraryInsertOrientationProbs{
  struct IOtriplet IOs[2]; // 1 : Oriented <- -> ....... 0 : Oriented -> <-; 2
};

struct placement{
  int assemPart;
  int offset1; // index on the assembly that this placement corresponds to
  int offset2;
  double likelihood; // likelihood of this placement (can be recalculated)
  unsigned char placeInfo; // orientation of this placement
  // placeInfo :
    // 1248
    // 00XXXXXX => Oriented <- ->
    // 10XXXXXX => Oriented -> <-
    // 11XXXXXX => Oriented -> ->
    // 01XXXXXX => Oriented <- <-
    // XX1XXXXX => This is the sequence on the assembly
    // XX0XXXXX => The compliment is the sequence on the assmebly
    // XXX1XXXX => left seq first
};

// the read_sequence data structure
// it contains the sequence and the sequence length
struct readSequence{
    unsigned char sequence[MAX_READ_LENGTH/4];
    char qval[MAX_READ_LENGTH];
    unsigned int seqLen;
};

// a paired read! takes (MAX_READ_LENGTH/2 + 5) bytes each.
struct pairedRead{
    struct readSequence readSubSeq[2]; // first sequence,second sequence
    float mu; // expected distance between sequences
    float sigma; // standard deviation of distances between sequences
    unsigned char readInfo; // info on stuff like orientation, qvals etc
    // readInfo :
    // 1248
    // 01XXXXXX => Oriented <- ->
    // 10XXXXXX => Oriented -> <-
    // 11XXXXXX => Oriented -> ->
    // 00XXXXXX => Oriented <- <-
    // XX1XXXXX => This is the sequence on the assembly
    // XX0XXXXX => The compliment is the sequence on the assmebly
    // XXX1XXXX => Qvals are present
    // XXX0XXXX => No Qvals
    // XXXX1XXX => Orientation is known
    // XXXX0XXX => Orientation is unknown
    // XXXXX1XX => compliment/non-compliment known/unknown
    unsigned char hasQval; // does the read have quality values?
    
    double placementNormalizer;
    
    struct placement placements[N_PLACEMENTS];
    unsigned char numPlacements;
};

struct assemblyPart{
  //char sequence[MAX_CONTIG_LENGTH];
  unsigned int seqLen;
  char *sequence;
  float *depthInward;
  float *depthOutward;
};

struct assembly{
  //struct assemblyPart assemblyParts[NUM_ASSEMBLY_PARTS];
  struct assemblyPart *assemblyParts;
};

typedef struct pairedRead pairedRead_t;
typedef struct readSequence readSequence_t;
typedef struct assemblyPart assemblyPart_t;
typedef struct assembly assembly_t;
typedef struct placement placement_t;
typedef struct libraryInsertOrientationProbs libInsOrProbs_t;

/***************************************************************
**** HELPER FUNCTIONS FOR WRAPPING/RE-WRAPPING THE SEQUENCE ****
***************************************************************/

unsigned char seqToChar(const char pos1, const char pos2, const char pos3, const char pos4);
void charToSeq(unsigned char num, char seq[], const int len);
pairedRead_t FindBestPlacements(pairedRead_t theRead, assemblyPart_t thePart, libInsOrProbs_t libProbs);
double GetInsertProbNormal(const double sigma, const double point);
double PlacementLikelihoodFTQ(const readSequence_t read, const assemblyPart_t part, const int offset);
double PlacementLikelihoodFCQ(const readSequence_t read, const assemblyPart_t part, const int offset);
double PlacementLikelihoodBTQ(const readSequence_t read, const assemblyPart_t part, const int offset);
double PlacementLikelihoodBCQ(const readSequence_t read, const assemblyPart_t part, const int offset);
double GetLikelihoodOfPlacementByErrors(const double expErrors);
double GetExpectedError(const char readRes, const char partRes, const double qval);
void combinePlacements(pairedRead_t theRead, placement_t tempPlacements[]);
int SuggestPlacement(placement_t oldPlacements[], placement_t newPlacement, unsigned char numPlacements);
void charToSeqFour(unsigned char num, char seq[]);
double lnfact(double input);

double getQtoP(char qualChar, int qOff) {
	int idx = qualChar - qOff;
	if (idx < 0 || idx >= 63 )
		printf("WARNING: getQtoP called out of range: %c %d %d\n", qualChar, qOff, idx);
	assert(idx >= 0 && idx < 63);
	return QtoP[idx];
}

void IncreaseAssemblyPartsByOne(assembly_t *theAssembly, int numParts){
  assemblyPart_t *tempPartPointer = malloc(numParts* sizeof(assemblyPart_t));
  int i;
  for(i = 0; i < numParts - 1; i++){
    tempPartPointer[i] = theAssembly->assemblyParts[i];
  }
  //free(theAssembly->assemblyParts);
  theAssembly->assemblyParts = tempPartPointer;
  //return theAssembly;
}

double poissonInt(int k, double lambda){
  int i;
  double answer = 1.0;
  for(i = 1; i < k+1; i++){
    answer = answer*lambda*exp(-lambda/(float)k)/(float)i;
  }
  return answer;
}

//uses Stirlings approximation to high precision
double lnfact(double input){
  return (input - 0.5)*log(input) - input + lnfactconst - 1.0/(12.0*input) - 1.0/(360.0*input*input*input) - 1.0/(1260.0*input*input*input*input*input);
}

// convert a 4-mer into a byte
unsigned char seqToChar(const char pos1, const char pos2, const char pos3, const char pos4){
    unsigned char seqChar = 0;
    seqChar += (pos1 == 'A' || pos1 == 'G');
    seqChar += (pos1 == 'A' || pos1 == 'C')*2;
    seqChar += (pos2 == 'A' || pos2 == 'G')*4;
    seqChar += (pos2 == 'A' || pos2 == 'C')*8;
    seqChar += (pos3 == 'A' || pos3 == 'G')*16;
    seqChar += (pos3 == 'A' || pos3 == 'C')*32;
    seqChar += (pos4 == 'A' || pos4 == 'G')*64;
    seqChar += (pos4 == 'A' || pos4 == 'C')*128;
    return seqChar;
}

// saves the quality information of a sequence
void makeQual(const char qual[], const unsigned int seqLen, char quality[], int qualityHistogram[]){
  int i;
  for(i = 0; i < seqLen; i++){
    quality[i] = qual[i];
    qualityHistogram[qual[i] - 64]++;
  }
}

// saves the quality information of a sequence pretending it is as high as it can be
void makeQualPerfect(const unsigned int seqLen, char quality[]){
  int i;
  for(i = 0; i < seqLen; i++){
    quality[i] = 127;
  }
}

// saves the quality information of a sequence
void makeAssemblySeq(const char seq[], const unsigned int seqLen, unsigned char sequence[]){
  int i;
  for(i = 0; i < seqLen; i++){
    sequence[i] = toupper(seq[i]);
  }
}

// converts a sequence to the condensed char equivalent
void makeSeq(const char seq[], const unsigned int seqLen, unsigned char sequence[]){
    int i;
    // make all multiples of 4 possible
    for(i = 0; i < seqLen/4; i++){
        sequence[i] = seqToChar(seq[i*4], seq[i*4 + 1], seq[i*4 + 2], seq[i*4 + 3]);
    }
    if(seqLen%4 == 3){
        sequence[seqLen/4] = seqToChar(seq[i*4], seq[i*4 + 1], seq[i*4 + 2], 'T');
    }else if(seqLen%4 == 2){
        sequence[seqLen/4] = seqToChar(seq[i*4], seq[i*4 + 1], 'T', 'T');
    }else if(seqLen%4 == 1){
        sequence[seqLen/4] = seqToChar(seq[i*4], 'T', 'T', 'T');
    }
}

// gives you the nucleotide at position loc from the sequence <seq>, not efficient for many close
// calls, use charToSeq to get 4 at a time.
char getCharFromSeqByLoc(const unsigned char seq[], const unsigned int loc){
    char subSeq[4];
    charToSeqFour(seq[loc/4], subSeq); // grab the 4-mer that we want
    return subSeq[loc%4]; // return the actual residue we want
}

static const char theFourConverter[256][4] = {"TTTT","GTTT","CTTT","ATTT","TGTT","GGTT","CGTT","AGTT","TCTT","GCTT","CCTT","ACTT","TATT","GATT","CATT","AATT","TTGT","GTGT","CTGT","ATGT","TGGT","GGGT","CGGT","AGGT","TCGT","GCGT","CCGT","ACGT","TAGT","GAGT","CAGT","AAGT","TTCT","GTCT","CTCT","ATCT","TGCT","GGCT","CGCT","AGCT","TCCT","GCCT","CCCT","ACCT","TACT","GACT","CACT","AACT","TTAT","GTAT","CTAT","ATAT","TGAT","GGAT","CGAT","AGAT","TCAT","GCAT","CCAT","ACAT","TAAT","GAAT","CAAT","AAAT","TTTG","GTTG","CTTG","ATTG","TGTG","GGTG","CGTG","AGTG","TCTG","GCTG","CCTG","ACTG","TATG","GATG","CATG","AATG","TTGG","GTGG","CTGG","ATGG","TGGG","GGGG","CGGG","AGGG","TCGG","GCGG","CCGG","ACGG","TAGG","GAGG","CAGG","AAGG","TTCG","GTCG","CTCG","ATCG","TGCG","GGCG","CGCG","AGCG","TCCG","GCCG","CCCG","ACCG","TACG","GACG","CACG","AACG","TTAG","GTAG","CTAG","ATAG","TGAG","GGAG","CGAG","AGAG","TCAG","GCAG","CCAG","ACAG","TAAG","GAAG","CAAG","AAAG","TTTC","GTTC","CTTC","ATTC","TGTC","GGTC","CGTC","AGTC","TCTC","GCTC","CCTC","ACTC","TATC","GATC","CATC","AATC","TTGC","GTGC","CTGC","ATGC","TGGC","GGGC","CGGC","AGGC","TCGC","GCGC","CCGC","ACGC","TAGC","GAGC","CAGC","AAGC","TTCC","GTCC","CTCC","ATCC","TGCC","GGCC","CGCC","AGCC","TCCC","GCCC","CCCC","ACCC","TACC","GACC","CACC","AACC","TTAC","GTAC","CTAC","ATAC","TGAC","GGAC","CGAC","AGAC","TCAC","GCAC","CCAC","ACAC","TAAC","GAAC","CAAC","AAAC","TTTA","GTTA","CTTA","ATTA","TGTA","GGTA","CGTA","AGTA","TCTA","GCTA","CCTA","ACTA","TATA","GATA","CATA","AATA","TTGA","GTGA","CTGA","ATGA","TGGA","GGGA","CGGA","AGGA","TCGA","GCGA","CCGA","ACGA","TAGA","GAGA","CAGA","AAGA","TTCA","GTCA","CTCA","ATCA","TGCA","GGCA","CGCA","AGCA","TCCA","GCCA","CCCA","ACCA","TACA","GACA","CACA","AACA","TTAA","GTAA","CTAA","ATAA","TGAA","GGAA","CGAA","AGAA","TCAA","GCAA","CCAA","ACAA","TAAA","GAAA","CAAA", "AAAA"};

void charToSeqFour(unsigned char num, char seq[]){
  seq[0] = theFourConverter[num][0];
  seq[1] = theFourConverter[num][1];
  seq[2] = theFourConverter[num][2];
  seq[3] = theFourConverter[num][3];
}

// unwraps a byte <num> into it's <len> residues and stores them in <seq>
// useful for output or direct manipulation of the sequence
void charToSeq(unsigned char num, char seq[], const int len){
    int i;
    for(i = 0; i < len; i++){
        if(num%2){ // A or G
            num = num >> 1;
            if(num%2){ // A or C match?
                seq[i] = 'A';
            }else{
                seq[i] = 'G';
            }
        }else{
            num = num >> 1;
            if(num%2){ // A or C match?
                seq[i] = 'C';
            }else{
                seq[i] = 'T';
            }
        }
        num = num >> 1;
    }
    for(i = len; i < 4; i++){ // fill in the rest with spaces
        seq[i] = ' ';
    }
}

// prints out the sequence in it's numeric (byte-wrapped) form
void PrintSequenceNumeric(const unsigned char sequence[], const unsigned int seqLen){
    int j;
    for(j = 0; j < seqLen/4; j++){
        printf("%i ", sequence[j]);
    }
    if(seqLen%4 != 0){ //print the end
        printf("%i", sequence[seqLen/4]);
    }
    printf("\n");
}

void PrintSequenceB(const unsigned char sequence[], const unsigned int seqLen){
    int j;
    for(j = seqLen-1; j > 0; j--){
        printf("%c", getCharFromSeqByLoc(sequence, j));
    }
    printf("\n");
}

void PrintSequence(const unsigned char sequence[], const unsigned int seqLen){
    char seqer[4];
    int j;
    for(j = 0; j < seqLen/4; j++){
        charToSeq(sequence[j], seqer, 4);
        printf("%.4s", seqer);
    }
    if(seqLen%4 != 0){
        charToSeq(sequence[seqLen/4], seqer, seqLen%4);
        printf("%.4s", seqer);
    }
    printf("\n");
}

void PrintQuality(const char quality[], const unsigned int seqLen){
  int i;
  for(i = 0; i < seqLen; i++){
    printf("%c", quality[i]);
  }
  printf("\n");
}

void PrintAssembly(const char sequence[], const unsigned int seqLen){
  int i;
  for(i = 0; i < seqLen; i++){
    printf("%c", sequence[i]);
  }
  printf("\n");
}

// get the numerical value for the quality of a base call
double getQualityP(const char quality[], const unsigned int i){
  return QtoP[quality[i] - 64];
}

/* This is a secret function, its magics are UNKNOWN */
int intMax(int a, int b){
    if(a > b){
        return a;
    }
    return b;
}

/* This is a secret function, its magics are UNKNOWN */
int intMin(int a, int b){
    if(a < b){
        return a;
    }
    return b;
}

char getComplimentRes(const char res){
  if(res == 'A'){return 'T';}
  if(res == 'T'){return 'A';}
  if(res == 'G'){return 'C';}
  return 'G';
}

void PrintPlacements(pairedRead_t theRead){
  int i;
  for(i = 0; i < theRead.numPlacements; i++){
    printf("L: %f, o1: %i, o2: %i, i: %i an: %i\n", theRead.placements[i].likelihood, theRead.placements[i].offset1, theRead.placements[i].offset2, theRead.placements[i].placeInfo, theRead.placements[i].assemPart);
  }
}

#endif
