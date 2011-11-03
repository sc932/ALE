// ALEhelpers.h

#ifndef _ALE_HELPERS_H_
#define _ALE_HELPERS_H_

#include <zlib.h>
#include <assert.h>
#include <sam.h>
#include <bam.h>
#include <search.h>
#include "kseq.h"
#include <math.h>

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

#define mapLens_MAX 25000
#define GCmaps_MAX 400
#define MAX_NAME_LENGTH 127

//static unsigned char powup[8] = {1,2,4,8,16,32,64,128};

// http://en.wikipedia.org/wiki/FASTQ_format
// generated from scripts.py
// lookup of q -> p. p == probability of being correct ( redefine Q=0 as 1/100, could argue for 0.25 for both Q=0 and 1 as there are only 4 choices )
static double QtoP[63] = {0.01,0.205671765276,0.36904265552,0.498812766373,0.601892829447,0.683772233983,0.748811356849,0.800473768503,0.841510680754,0.874107458821,0.9,0.920567176528,0.936904265552,0.949881276637,0.960189282945,0.968377223398,0.974881135685,0.98004737685,0.984151068075,0.987410745882,0.99,0.992056717653,0.993690426555,0.994988127664,0.996018928294,0.99683772234,0.997488113568,0.998004737685,0.998415106808,0.998741074588,0.999,0.999205671765,0.999369042656,0.999498812766,0.999601892829,0.999683772234,0.999748811357,0.999800473769,0.999841510681,0.999874107459,0.9999,0.999920567177,0.999936904266,0.999949881277,0.999960189283,0.999968377223,0.999974881136,0.999980047377,0.999984151068,0.999987410746,0.99999,0.999992056718,0.999993690427,0.999994988128,0.999996018928,0.999996837722,0.999997488114,0.999998004738,0.999998415107,0.999998741075,0.999999,0.999999205672,0.999999369043};
// lookup of q -> log(p)
static double QtoLogP[63] = {-4.60517018598809,-1.58147375340709,-0.996843044007323,-0.695524471331768,-0.507675873695919,-0.380130408066409,-0.289268187201664,-0.222551515972973,-0.172556572913703,-0.134551960288297,-0.105360515657826,-0.0827653026686952,-0.0651741731994102,-0.0514182741582494,-0.0406248442212775,-0.0321335740234123,-0.0254397275341249,-0.0201543647615499,-0.0159758692470951,-0.0126691702086946,-0.0100503358535015,-0.00797499827826794,-0.00632956293131057,-0.00502447389071099,-0.00398901726687351,-0.0031672882259884,-0.00251504651167362,-0.0019972555025822,-0.00158615046381827,-0.00125971852431254,-0.00100050033358353,-0.000794643880834491,-0.000631156481354146,-0.000501312870301668,-0.000398186436697967,-0.000316277776543464,-0.000251220196151159,-0.000199546139006532,-0.000158501879759277,-0.000125900466131129,-0.000100005000333347,-7.94359779538046e-05,-6.30977246195429e-05,-5.01199789851444e-05,-3.9811509467615e-05,-3.16232770105271e-05,-2.51191794839762e-05,-1.99528220562603e-05,-1.58490575956675e-05,-1.25893332453737e-05,-1.00000500002878e-05,-7.94331354807621e-06,-6.30959290540238e-06,-5.01188455944388e-06,-3.98107992443596e-06,-3.16228300002451e-06,-2.51188915476935e-06,-1.99526399049202e-06,-1.58489425589492e-06,-1.25892579247174e-06,-1.00000050002909e-06,-7.94328315520852e-07,-6.30957199034177e-07};
// lookup of q -> log((1.0 - p)/3.0)
// assume all remaining probability goes evenly to remainder of possible choices
static double QtoLogPMiss[63] = {-1.10866262452161,-1.32887079796787,-1.55912930726722,-1.78938781656687,-2.01964632586698,-2.24990483516462,-2.48016334446437,-2.71042185376338,-2.94068036306405,-3.17093887236606,-3.40119738166216,-3.63145589096695,-3.86171440026127,-4.09197290955493,-4.32223141886856,-4.55248992814918,-4.7827484374624,-5.01300694674239,-5.24326545603285,-5.47352396535216,-5.7037824746562,-5.93404098398617,-6.16429949322361,-6.39455800260883,-6.624816511737,-6.85507502120647,-7.08533353025741,-7.31559203973642,-7.54585054934241,-7.77610905818733,-8.00636756765025,-8.23662607660251,-8.46688458701014,-8.69714309480488,-8.92740160372637,-9.15766011420038,-9.38791862404765,-9.61817713523669,-9.84843564359849,-10.0786941527697,-10.3089526606444,-10.5392111758906,-10.769469686344,-10.9997281957802,-11.2299866992325,-11.4602451945463,-11.6905037289838,-11.9207622332408,-12.1510207302817,-12.3812792537034,-12.6115377536429,-12.8417963066474,-13.0720548427394,-13.3023133486375,-13.5325717168828,-13.7628301926671,-13.9930889812304,-14.2233474766155,-14.4536059494995,-14.6838646644129,-14.9141228466036,-15.1443816513791,-15.3746404112262};

const static double lnfactconst = 0.918938533204672741780329;

static const char WELCOME_MSG[80] = "Welcome to the Assembly Likelihood Estimator!\n(C) 2010 Scott Clark\n\n";
static const char USAGE[80] = "Usage: %s [-options] readSorted.[s|b]am assembly.fasta[.gz] ALEoutput.txt\n";
static const char SHORT_OPTIONS[80] = "    Options:\n    -h : print out help\n";
static const char LONG_OPTIONS[1024] = "Options: <i>nt <f>loat [default]\n  -h        : print out this help\n\n   -kmer <f> : Kmer depth for kmer stats [4]\n  -qOff <i> : Quality ascii offset (illumina) [33] or 64 (or 0)\n  -pl placementOutputBAM\n\n\n";

static const char theFourConverter[256][4] = {"TTTT","GTTT","CTTT","ATTT","TGTT","GGTT","CGTT","AGTT","TCTT","GCTT","CCTT","ACTT","TATT","GATT","CATT","AATT","TTGT","GTGT","CTGT","ATGT","TGGT","GGGT","CGGT","AGGT","TCGT","GCGT","CCGT","ACGT","TAGT","GAGT","CAGT","AAGT","TTCT","GTCT","CTCT","ATCT","TGCT","GGCT","CGCT","AGCT","TCCT","GCCT","CCCT","ACCT","TACT","GACT","CACT","AACT","TTAT","GTAT","CTAT","ATAT","TGAT","GGAT","CGAT","AGAT","TCAT","GCAT","CCAT","ACAT","TAAT","GAAT","CAAT","AAAT","TTTG","GTTG","CTTG","ATTG","TGTG","GGTG","CGTG","AGTG","TCTG","GCTG","CCTG","ACTG","TATG","GATG","CATG","AATG","TTGG","GTGG","CTGG","ATGG","TGGG","GGGG","CGGG","AGGG","TCGG","GCGG","CCGG","ACGG","TAGG","GAGG","CAGG","AAGG","TTCG","GTCG","CTCG","ATCG","TGCG","GGCG","CGCG","AGCG","TCCG","GCCG","CCCG","ACCG","TACG","GACG","CACG","AACG","TTAG","GTAG","CTAG","ATAG","TGAG","GGAG","CGAG","AGAG","TCAG","GCAG","CCAG","ACAG","TAAG","GAAG","CAAG","AAAG","TTTC","GTTC","CTTC","ATTC","TGTC","GGTC","CGTC","AGTC","TCTC","GCTC","CCTC","ACTC","TATC","GATC","CATC","AATC","TTGC","GTGC","CTGC","ATGC","TGGC","GGGC","CGGC","AGGC","TCGC","GCGC","CCGC","ACGC","TAGC","GAGC","CAGC","AAGC","TTCC","GTCC","CTCC","ATCC","TGCC","GGCC","CGCC","AGCC","TCCC","GCCC","CCCC","ACCC","TACC","GACC","CACC","AACC","TTAC","GTAC","CTAC","ATAC","TGAC","GGAC","CGAC","AGAC","TCAC","GCAC","CCAC","ACAC","TAAC","GAAC","CAAC","AAAC","TTTA","GTTA","CTTA","ATTA","TGTA","GGTA","CGTA","AGTA","TCTA","GCTA","CCTA","ACTA","TATA","GATA","CATA","AATA","TTGA","GTGA","CTGA","ATGA","TGGA","GGGA","CGGA","AGGA","TCGA","GCGA","CCGA","ACGA","TAGA","GAGA","CAGA","AAGA","TTCA","GTCA","CTCA","ATCA","TGCA","GGCA","CGCA","AGCA","TCCA","GCCA","CCCA","ACCA","TACA","GACA","CACA","AACA","TTAA","GTAA","CTAA","ATAA","TGAA","GGAA","CGAA","AGAA","TCAA","GCAA","CCAA","ACAA","TAAA","GAAA","CAAA", "AAAA"};

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

struct contig_struct{
  char *name;
  int seqLen;
  unsigned char isCircular;
  char *seq;
  float *depth;
  float *matchLikelihood;
  float *depthLikelihood;
  float *kmerLikelihood;
  unsigned char *GCcont; // range of 0 - 100
};

struct assembly_struct{
  int numContigs;
  long totalAssemLen;
  struct contig_struct **contigs;
};

struct setOfAlignments{
  double likelihood;
  int start1, start2;
  int end1, end2;
  int contigId;
  char name[MAX_NAME_LENGTH+1];
  struct setOfAlignments *nextAlignment;
};

enum MATE_ORIENTATION {
  VALID_FR,
  VALID_RF,
  VALID_FF,
  NOT_PROPER_FR,
  NOT_PROPER_RF,
  NOT_PROPER_FF,
  CHIMER,
  READ1_ONLY, // but is paired
  READ2_ONLY, // but is paired
  HALF_VALID_MATE, // paired but only one read is observed
  SINGLE_READ,// not paired
  NO_READS,
  UNRELATED_PAIR, // two reads, each paired, but not to each other (i.e. not in sort-by-name order)
  UNMAPPED_SINGLE, // not paired, not mapped
  UNMAPPED_PAIR,   // paired, neither mapped
  MATE_ORIENTATION_MAX
};

const static char *MATE_ORIENTATION_LABELS[MATE_ORIENTATION_MAX] = {
  "FR",
  "RF",
  "FF",
  "NOT_PROPER_FR",
  "NOT_PROPER_RF",
  "NOT_PROPER_FF",
  "CHIMER",
  "READ1_ONLY",
  "READ2_ONLY",
  "HALF_VALID_MATE",
  "SINGLE_READ",
  "NO_READS",
  "UNRELATED_PAIR",
  "UNMAPPED_SINGLE",
  "UNMAPPED_PAIR"
};

struct libraryMateParameters {
  double insertLength; // insert length mean
  double insertStd; // insert length std dev
  double libraryFraction; // fraction of mates that map in this orientation
  long count;
  int isValid;
};

typedef struct libraryMateParameters libraryMateParametersT;

struct libraryParameters {
  libraryMateParametersT mateParameters[MATE_ORIENTATION_MAX];
  long avgReadSize;
  long numReads;
  double totalValidSingleFraction;
  double totalValidMateFraction;
  double totalChimerMateFraction;
  double totalUnmappedFraction;
  int qOff;
  int isSortedByName;
};

typedef struct setOfAlignments alignSet_t;
typedef struct contig_struct contig_t;
typedef struct assembly_struct assemblyT;
typedef struct libraryParameters libraryParametersT;

struct _contig_ll {
  contig_t *contig;
  void *next;
};
typedef struct _contig_ll contig_ll;


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
double getQtoP(char qualChar, int qOff);
double getQtoLogP(char qualChar, int qOff);
double getQtoLogPMiss(char qualChar, int qOff);
void IncreaseAssemblyPartsByOne(assembly_t *theAssembly, int numParts);
double poissonInt(int k, double lambda);

void initAlignment(alignSet_t *dst);

void copyAlignment(alignSet_t *dst, const alignSet_t *src);

void swap(void **x, void **y);

void printAssembly(assemblyT *theAssembly);

int isGC(char seq);

int getGCtotal(char seq1[], int seq1len);

int findNumAssemPieces(kseq_t *ins);

void readAssembly(kseq_t *ins, assemblyT *theAssembly);

// below is my attempt at a hanning window convolution, I coded it from scratch so watch for bugs!
void calculateGCcont(assemblyT *theAssembly, int windowSize);

int getSeqLenBAM(bam1_t *read);

int getMapLenBAM(bam1_t *read1);

enum MATE_ORIENTATION getPairedMateOrientation(bam1_t *read1);

enum MATE_ORIENTATION getMateOrientation(bam1_t *read1, bam1_t *read2);

enum MATE_ORIENTATION readMatesBAM(samfile_t *ins, libraryParametersT *libParams, bam1_t *read1, bam1_t *read2);

int assemblySanityCheck(assemblyT *theAssembly);

assemblyT *loadAssembly(char *filename);

void freeContig(contig_t *contig);

void freeAssembly(assemblyT *theAssembly);

samfile_t *openSamOrBam(const char *fileName);
#endif
