// ALEhelpers.c

#include "ALEhelpers.h"

double getQtoP(char qualChar, int qOff) {
    int idx = qualChar - qOff;
    if (idx < 0 || idx >= 63 )
        //printf("WARNING: getQtoP called out of range: %c %d %d\n", qualChar, qOff, idx);
    assert(idx >= 0 && idx < 63);
    return QtoP[idx];
}
double getQtoLogP(char qualChar, int qOff) {
    int idx = qualChar - qOff;
    if (idx < 0 || idx >= 63 )
        //printf("WARNING: getQtoLogP called out of range: %c %d %d\n", qualChar, qOff, idx);
    assert(idx >= 0 && idx < 63);
    return QtoLogP[idx];
}
double getQtoLogPMiss(char qualChar, int qOff) {
    int idx = qualChar - qOff;
    if (idx < 0 || idx >= 63 )
        //printf("WARNING: getQtoLogPMiss called out of range: %c %d %d\n", qualChar, qOff, idx);
    assert(idx >= 0 && idx < 63);
    return QtoLogPMiss[idx]; // TODO switch to (1-Q)*Q
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
        //printf("%i ", sequence[j]);
    }
    if(seqLen%4 != 0){ //print the end
        //printf("%i", sequence[seqLen/4]);
    }
    //printf("\n");
}

void PrintSequenceB(const unsigned char sequence[], const unsigned int seqLen){
    int j;
    for(j = seqLen-1; j > 0; j--){
        //printf("%c", getCharFromSeqByLoc(sequence, j));
    }
    //printf("\n");
}

void PrintSequence(const unsigned char sequence[], const unsigned int seqLen){
    char seqer[4];
    int j;
    for(j = 0; j < seqLen/4; j++){
        charToSeq(sequence[j], seqer, 4);
        //printf("%.4s", seqer);
    }
    if(seqLen%4 != 0){
        charToSeq(sequence[seqLen/4], seqer, seqLen%4);
        //printf("%.4s", seqer);
    }
    //printf("\n");
}

void PrintQuality(const char quality[], const unsigned int seqLen){
  int i;
  for(i = 0; i < seqLen; i++){
    //printf("%c", quality[i]);
  }
  //printf("\n");
}

void PrintAssembly(const char sequence[], const unsigned int seqLen){
  int i;
  for(i = 0; i < seqLen; i++){
    //printf("%c", sequence[i]);
  }
  //printf("\n");
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
    //printf("L: %f, o1: %i, o2: %i, i: %i an: %i\n", theRead.placements[i].likelihood, theRead.placements[i].offset1, theRead.placements[i].offset2, theRead.placements[i].placeInfo, theRead.placements[i].assemPart);
  }
}

void initAlignment(alignSet_t *dst) {
    dst->likelihood = 0.0;
    dst->start1 = -1;
    dst->start2 = -1;
    dst->end1 = -1;
    dst->end2 = -1;
    dst->contigId1 = -1;
    dst->contigId2 = -1;
    dst->name = NULL;
    dst->nextAlignment = NULL;
}

void destroyAlignment(alignSet_t *dst) {
    if (dst->name != NULL)
        free(dst->name);
}

void copyAlignment(alignSet_t *dst, const alignSet_t *src) {
    assert(dst != NULL);
    assert(src != NULL);
    destroyAlignment(dst);

    dst->likelihood = src->likelihood;
    dst->start1 = src->start1;
    dst->start2 = src->start2;
    dst->end1 = src->end1;
    dst->end2 = src->end2;
    dst->contigId1 = src->contigId1;
    dst->contigId2 = src->contigId2;
    dst->name = strdup(src->name);
    dst->nextAlignment = src->nextAlignment;
}

void swap(void **x, void **y) {
    void *t = *x;
    *x = *y;
    *y = t;
}

int isGC(char seq) {
    return (seq == 'G' || seq == 'C' || seq == 'g' || seq == 'c');
}

int getGCtotal(char seq1[], int seq1len){
    int GCtot = 0, i;
    for(i = 0; i < seq1len; i++){
        if(isGC(seq1[i])){
            GCtot++;
        }
    }
    return GCtot;
}

int findNumAssemPieces(kseq_t *ins){
    int l, count = 0;
    while ((l = kseq_read(ins)) >= 0) {
        count++;
    }
    return count;
}

void readAssembly(kseq_t *ins, assemblyT *theAssembly){
    int contigLen, i, l, j = 0;

    // read contigs into a linked list.  Consolidate to an array later
    contig_ll *tmp, *head = malloc(sizeof(contig_ll));
    head->contig = NULL;
    head->next = NULL;
    tmp = head;

    while ((l = kseq_read(ins)) >= 0) {
        contigLen = (int)(ins->seq.l);
        ////printf("Found contig %d, contigLen = %i, name=%s\n", j, contigLen, ins->name.s);

        contig_t *contig = tmp->contig = (contig_t*) malloc(sizeof(contig_t));
        contig->seqLen = contigLen;
        contig->isCircular = 0; // assume nothing is circular until there reasonable evidence
        contig->seq = malloc(contigLen*sizeof(char));
        contig->depth = malloc(contigLen*sizeof(float));
        contig->matchLikelihood = malloc(contigLen*sizeof(float));
        contig->depthLikelihood = malloc(contigLen*sizeof(float));
        contig->kmerLikelihood = malloc(contigLen*sizeof(float));
        contig->GCcont = malloc(contigLen*sizeof(unsigned char));
        contig->name = strdup(ins->name.s);
        for(i = 0; i < contigLen; i++){
            contig->seq[i] = toupper(ins->seq.s[i]);
            contig->depth[i] = 0.0;
            contig->matchLikelihood[i] = 0.0;
            contig->depthLikelihood[i] = 0.0;
            contig->kmerLikelihood[i] = 0.0;
            contig->GCcont[i] = 0;
        }
        j++;
        tmp->next = malloc(sizeof(contig_ll));
        tmp = tmp->next;
        tmp->contig = NULL;
        tmp->next = NULL;

    }
    int numberAssemblyPieces = j;
    //printf("Found %d contigs\n", numberAssemblyPieces);

    theAssembly->contigs = (contig_t**)malloc((numberAssemblyPieces)*sizeof(contig_t*));
    theAssembly->numContigs = numberAssemblyPieces;
    theAssembly->totalScore = 0.0;
    theAssembly->kmerAvgSum = 0.0;
    theAssembly->kmerAvgNorm = 0.0;
    theAssembly->depthScoreAvgSum = 0.0;
    theAssembly->depthScoreAvgNorm = 0.0;
    theAssembly->depthAvgSum = 0.0;
    theAssembly->depthAvgNorm = 0.0;
    theAssembly->totalUnmappedReads = 0;

    // consolidate linked list into array, free linked list
    tmp = head;
    i = 0;
    while (tmp != NULL) {
        if (tmp->contig != NULL)
            theAssembly->contigs[i++] = tmp->contig;
        head = tmp;
        tmp = head->next;
        free(head);
    }
}

int validateAssemblyIsSameAsAlignment(bam_header_t *header, assemblyT *theAssembly) {
    int i;
    //printf("Validating assembly and alignment files consisting of %d contigs\n", header->n_targets);
    if (header->n_targets != theAssembly->numContigs) {
        //printf("Different number of contigs in assembly (%d) and alignmentfile (%d)\n", theAssembly->numContigs, header->n_targets);
        return 0;
    }
    for (i = 0; i < header->n_targets; i++) {
        if (header->target_len[i] != theAssembly->contigs[i]->seqLen) {
        	//printf("Different contig length in contig %d: %d vs %d\n", i, header->target_len[i], theAssembly->contigs[i]->seqLen);
        	return 0;
        }
        if (strcmp(header->target_name[i], theAssembly->contigs[i]->name) != 0) {
        	//printf("Warning assembly and alignment files disagree on the name of contig %d: %s vs %s\n", i, header->target_name[i], theAssembly->contigs[i]->name);
        }
    }
    return 1;
}

// below is my attempt at a hanning window convolution, I coded it from scratch so watch for bugs!
void calculateGCcont(assemblyT *theAssembly, int windowSize){
    int i, j, baseGC;
    int *GCpast = malloc(sizeof(int) * windowSize);
    int skippedContigs = 0;
    for(i = 0; i < theAssembly->numContigs; i++){
        contig_t *contig = theAssembly->contigs[i];
        if (contig->seqLen < 2 * windowSize) {
        	for(j=0; j < contig->seqLen ; j++) {
        		contig->GCcont[j] = 101; // throw away
        	}
        	skippedContigs++;
        	continue; // contig is too small to process
        }
        baseGC = getGCtotal(contig->seq, windowSize);
        GCpast[0] = baseGC;
        for(j = 0; j < windowSize; j++){
        	contig->GCcont[j] = floor(100.0*(double)baseGC/(double)((j+1)*windowSize));
        	GCpast[(j+1)%windowSize] = GCpast[j%windowSize];
        	if(isGC(contig->seq[j])){
        		GCpast[(j+1)%windowSize]--;
        	}
        	if(isGC(contig->seq[j+windowSize])){
        		GCpast[(j+1)%windowSize]++;
        	}
        	baseGC += GCpast[(j+1)%windowSize];
        }
        for(j = windowSize; j < contig->seqLen - windowSize; j++){
        	contig->GCcont[j] = floor(100.0*(double)baseGC/(double)(windowSize*windowSize));
        	baseGC -= GCpast[(j+1)%windowSize];
        	GCpast[(j+1)%windowSize] = GCpast[j%windowSize];
        	if(isGC(contig->seq[j])){
        		GCpast[(j+1)%windowSize]--;
        	}
        	if(isGC(contig->seq[j+windowSize])) {
        		GCpast[(j+1)%windowSize]++;
        	}
        	baseGC += GCpast[(j+1)%windowSize];
        }
        for(j = contig->seqLen - windowSize; j < contig->seqLen; j++){
        	contig->GCcont[j] = floor(100.0*(double)baseGC/(double)((contig->seqLen - j)*windowSize));
        	baseGC -= GCpast[(j+1)%windowSize];
        }
    }
    free(GCpast);
    //printf("%d contigs were too small to calculate GC coverage over a %d window\n", skippedContigs, windowSize);
}

int getSeqMapLenBAM(bam1_t *read) {
    assert(read != NULL);
    return bam_cigar2qlen(&read->core, bam1_cigar(read));
}

int getFragmentMapLenBAM(bam1_t *read1) {
    assert(read1 != NULL);

    int left = read1->core.pos < read1->core.mpos ? read1->core.pos : read1->core.mpos;
    int readLength = getSeqMapLenBAM(read1);
    int right1 = read1->core.pos + readLength;
    int right2 = read1->core.mpos + readLength;
    int right = right1 > right2 ? right1 : right2;
    assert(right >= left);
    return right - left;
}

enum MATE_ORIENTATION getPairedMateOrientation(bam1_t *read1) {
    if ((read1->core.flag & BAM_FUNMAP) == BAM_FUNMAP || (read1->core.flag & BAM_FMUNMAP) == BAM_FMUNMAP) {
        // read or mate is not mapped
        if ((read1->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
        	// paired
        	if ((read1->core.flag & BAM_FUNMAP) == BAM_FUNMAP && (read1->core.flag & BAM_FMUNMAP) == BAM_FMUNMAP) {
        		// neither read is mapped
        		return UNMAPPED_PAIR;
        	} else {
        		// only one read in the pair is mapped
        		if ((read1->core.flag & BAM_FREAD1) == BAM_FREAD1) {
        			// this is READ1
        			return (read1->core.flag & BAM_FUNMAP) == BAM_FUNMAP ? READ2_ONLY : READ1_ONLY;
        		} else if ((read1->core.flag & BAM_FREAD2) == BAM_FREAD2) {
        			// this is READ2
        			return (read1->core.flag & BAM_FUNMAP) == BAM_FUNMAP ? READ1_ONLY : READ2_ONLY;
        		} else {
        			// not simply a pair of two reads, and this is not mapped
        			return UNMAPPED_SINGLE;
        		}
        	}
        } else {
        	// not mapped and not paired
        	return UNMAPPED_SINGLE;
        }
    }

    // this read is mapped
    assert((read1->core.flag & BAM_FUNMAP) != BAM_FUNMAP);

    if ((read1->core.flag & BAM_FPAIRED) != BAM_FPAIRED) {
        return SINGLE_READ;
    }
    if (((read1->core.flag & (BAM_FREAD1 | BAM_FREAD2)) == 0) || (read1->core.isize == 0)) {
        // required for PacBio reads that claim they are pairs but are not in a strict sense
        return SINGLE_READ;
    }
    assert((read1->core.flag & BAM_FPAIRED) == BAM_FPAIRED);

    char isProper = ((read1->core.flag & BAM_FPROPER_PAIR) == BAM_FPROPER_PAIR);
    if (read1->core.tid == read1->core.mtid) {
        // reads map to same contig
        int read1Dir = (read1->core.flag & BAM_FREVERSE) == BAM_FREVERSE ? 1 : 0;
        int read2Dir = (read1->core.flag & BAM_FMREVERSE) == BAM_FMREVERSE ? 1 : 0;
        if (read1Dir == read2Dir)
        	return isProper ? VALID_FF : NOT_PROPER_FF;
        else {
        	// TODO rethink this if read sizes are different could use read1->core.isize instead
        	int readLength = getSeqMapLenBAM(read1);
        	if (read1Dir == 0) {
        		if (read1->core.pos <= read1->core.mpos + readLength)
        			return isProper ? VALID_FR : NOT_PROPER_FR;
                else
        		    return isProper ? VALID_RF : NOT_PROPER_RF;
        	} else {
        		if (read1->core.mpos <= read1->core.pos + readLength)
        			return isProper ? VALID_FR : NOT_PROPER_FR;
                else
        		    return isProper ? VALID_RF : NOT_PROPER_RF;
        	}
        }
    } else {
        // reads map to different contig
        return CHIMER;
    }

}

enum MATE_ORIENTATION readNextBAM(samfile_t *ins, libraryParametersT *libParams, bam1_t *read1) {
    assert(ins != NULL);
    assert(read1 != NULL);

    int bytesRead = samread(ins, read1);
    if (bytesRead <= 0)
        return NO_READS;
    else
        return getPairedMateOrientation(read1);
}


// prints out all of the alignments in the linked list
void printAlignments(alignSet_t *head){
    // print out the head
    //printf("Alignment 1 for read %s: %f at %i-%i and %i-%i.\n", head->name, head->likelihood, head->start1, head->end1, head->start2, head->end2);
    alignSet_t *current = head;
    int i = 1;
    while(current->nextAlignment != NULL){
        current = current->nextAlignment;
        i++;
        //printf("Alignment %i for read %s: %f at %i-%i and %i-%i.\n", i, current->name, current->likelihood, current->start1, current->end1, current->start2, current->end2);
    }
}

void writeToOutput(assemblyT *theAssembly, FILE *out){
    int i, j;
    //printf("Writing statistics to output file.\n");
    fprintf(out, "# ALE_score: %lf\n", theAssembly->totalScore);
    fprintf(out, "# numContigs: %d\n# totalAssemLen: %ld\n# kmerAvg: %lf\n# depthScoreAvg: %lf\n# depthAvg: %lf\n# totalUnmappedReads: %d\n# readAvgLen: %lf\n", theAssembly->numContigs, theAssembly->totalAssemLen, theAssembly->kmerAvgSum/theAssembly->kmerAvgNorm, theAssembly->depthScoreAvgSum/theAssembly->depthScoreAvgNorm, theAssembly->depthAvgSum/theAssembly->depthAvgNorm, theAssembly->totalUnmappedReads, theAssembly->readAvgLen);
    // TODO TURNED OFF FOR ASSEMBLATHON
    /*
    for(i = 0; i < theAssembly->numContigs; i++){
        contig_t *contig = theAssembly->contigs[i];
        fprintf(out, "# Reference: %s %i\n# contig position depth ln(depthLike) ln(placeLike) ln(kmerLike) ln(totalLike)\n", contig->name, contig->seqLen);
        for(j = 0; j < contig->seqLen; j++){
            float logTotal = contig->depthLikelihood[j] + contig->matchLikelihood[j] + contig->kmerLikelihood[j];
            fprintf(out, "%d %d %0.3f %0.3f %0.3f %0.3f %0.3f\n", i, j, contig->depth[j], contig->depthLikelihood[j], contig->matchLikelihood[j], contig->kmerLikelihood[j], logTotal);
        }
    }
    */
    //printf("Total ALE Score: %lf\n", theAssembly->totalScore);
}

int assemblySanityCheck(assemblyT *theAssembly){
    int i, j, num = theAssembly->numContigs;
    int error = 1;
    for(j=0; j < num ; j++){
        contig_t *contig = theAssembly->contigs[j];
        for(i = 0; i < contig->seqLen; i++){
            if(contig->seq[i] != 'A' && contig->seq[i] != 'T' && contig->seq[i] != 'C' && contig->seq[i] != 'G' && contig->seq[i] != 'N'){
                printf("Found an error in the assembly, contig %d, position %d = %c\n", j, i, contig->seq[i]);
                contig->seq[i] = 'N';
                error = 0;
            }
        }
    }
    if(error == 0){
      printf("ALE considers these errors (%d) and will treat them as such; leaving a low depth, kmer score and placement likelihood around these regions. ALE only accepts the bases A,T,C,G,N.\n", error);
    }
    return error;
}

assemblyT *loadAssembly(char *filename) {

    // attempt to open the input file
    gzFile *assemblyFile = gzopen(filename, "r");
    kseq_t *Aseq;
    if(assemblyFile == NULL){
        //printf("Error! Could not open assembly file: %s\n", filename);
        exit(1);
    }

    assemblyT *theAssembly = malloc(sizeof(assemblyT));
    if (theAssembly == NULL)
        exit(1);

    Aseq = kseq_init(assemblyFile);

    readAssembly(Aseq, theAssembly);

    kseq_destroy(Aseq);

    //printf("Done reading in assembly.\n");

    //printAssembly(theAssembly);
    assemblySanityCheck(theAssembly);
    gzclose(assemblyFile);

    return theAssembly;
}

void freeContig(contig_t *contig) {
    if (contig == NULL)
        return;
    free(contig->name);
    free(contig->seq);
    free(contig->depth);
    free(contig->matchLikelihood);
    free(contig->kmerLikelihood);
    free(contig->depthLikelihood);
    free(contig->GCcont);
    free(contig);
}

void freeAssembly(assemblyT *theAssembly) {
    if (theAssembly != NULL) {
        if (theAssembly->contigs != NULL) {
        	int i;
        	for (i = 0; i < theAssembly->numContigs; i++)
        		freeContig(theAssembly->contigs[i]);
        	free(theAssembly->contigs);
        }
        free(theAssembly);
    }
}

samfile_t *openSamOrBam(const char *fileName) {
    samfile_t *in = samopen(fileName, "rb", 0);
    if (in == NULL || in->header == NULL) {
        //printf("Checking if %s is a SAM formatted file\n", fileName);
        in = samopen(fileName, "r", 0);
        if (in == NULL || in->header == NULL) {
            //printf("Error! Failed to open BAM/SAM file %s\n", fileName);
            exit(1);
        }
    }
    return in;
}

