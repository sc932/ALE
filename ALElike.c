// ALElike.c

#include "ALElike.h"

// casts a single numeric char to its int
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
  return 0;
}

//uses Stirlings approximation to high precision
double lnfact2(double input){
  return (input - 0.5)*log(input) - input + lnfactconst2 - 1.0/(12.0*input) - 1.0/(360.0*input*input*input) - 1.0/(1260.0*input*input*input*input*input);
}

// finds the log poisson pmf at value k for mean lambda
double poissonPMF(double k, double lambda){
  return k*log(lambda) - lambda - lnfact2(k + 1);
}

// finds the insert probability (assuming normal distribution) P(point | N(0,sigma))
double GetInsertProbNormal(double point, const double sigma){
  double p1 = erf((point + 0.5)/sqrt(2*sigma*sigma));
  double p2 = erf((point - 0.5)/sqrt(2*sigma*sigma));
  double prob = 0.5*(p1 - p2);
  //printf("Point: %lf, p1: %lf, p2: %lf = %lf\n", point, erf((point + 0.5)/sqrt(2*sigma*sigma)), erf((point - 0.5)/sqrt(2*sigma*sigma)), prob);
  return prob;
}

double getInsertLikelihoodBAM(bam1_t *read1, double mu, double var){
  assert(read1 != NULL);
  double mapLength = getMapLenBAM(read1);
  assert(mapLength > 0.0);
  double likelihood = GetInsertProbNormal(abs(mapLength - mu), var);
  //printf("getInsertLikelihoodBAM(%f,%f): %e\n", mu, var,likelihood);
  assert(likelihood >= 0.0);
  return likelihood;
}

// finds the loglikelihood of a string of misses in read from seqPos to matchLen
double loglikeMiss(char *readQual, int seqPos, int missLen, int qOff){
  int i;
  double loglikelihood = 0.0;
  for(i = seqPos; i < seqPos + missLen; i++){
    //likelihood = likelihood*((1.0 - 0.99)/3.0);
    //likelihood = likelihood*((1.0 - getQtoP(readQual[i], qOff))/3.0); // sometimes want -33/64?
	double logp = getQtoLogPMiss(readQual[i], qOff);
	//printf("likeMiss(%d, %d, %c): %lf %lf %lf %lf\n", missLen, i, readQual[i] + 33, logp, loglikelihood, exp(logp), exp(loglikelihood));
	loglikelihood += logp;
  }
  //printf("loglikeMiss(%d): %e\n", missLen, loglikelihood);
  assert(loglikelihood <= 0.0);
  return loglikelihood;
}

// finds the loglikelihood of a string of matches in read from seqPos to matchLen
double loglikeMatch(char *readQual, int seqPos, int matchLen, int qOff){
  int i;
  double loglikelihood = 0.0;
  for(i = seqPos; i < seqPos + matchLen; i++){
    loglikelihood += getQtoLogP(readQual[i], qOff);// sometimes want -33/64?
    //printf("likeMatch %d %d %f %f %f\n", i, readQual[i], QtoLogP[readQual[i] - qOff], getQtoLogP(readQual[i], qOff), loglikelihood);
  }
  //printf("loglikeMatch(%d): %e\n", matchLen, loglikelihood);
  assert(loglikelihood <= 0.0);
  return loglikelihood;
}

// finds the loglikelihood of an insertion (right now it is the same as a miss)
double loglikeInsertion(char *readQual, int seqPos, int insertionLength, int qOff) {
  // assume as unlikely as a substitution
  // TODO refine
  double loglikelihood = loglikeMiss(readQual, seqPos, insertionLength, qOff);
  //printf("loglikeInsertion(%d): %e\n", insertionLength, loglikelihood);
  return loglikelihood;
}

// finds the loglikelihood of an deletion (right now it is the same as a miss)
double loglikeDeletion(char *readQual, int seqPos, int deletionLength, int qOff) {
  // assume as unlikely as a substitution of previous base
  // TODO refine
  int delPos = (seqPos > 0) ? seqPos - 1 : seqPos;
  assert(delPos >= 0);
  double loglikelihood = loglikeMiss(readQual, delPos, 1, qOff) * (double)deletionLength;
  //printf("loglikeDeletion(%d): %e\n", deletionLength, loglikelihood);
  return loglikelihood;
}

// used to reduce loglikelihood in case of missmatches only
// (CIGAR already has accounted for matchlength, inserts, deletions)
double getMDLogLikelihood(char *MD, char *readQual, int qOff) {
  assert(MD != NULL && MD[0] != '\0');
  assert(readQual != NULL);
  //assert(readQual[0] != '\0');

  int stop = 0;
  int pos = 0;
  int seqPos = 0;
  double loglikelihood = 0.0;

  // parse MD field
  while(stop == 0){
    // matches
	int seqCount = 0;
	while(isdigit(MD[pos])){
	  	seqCount = seqCount*10 + (int)(MD[pos]) - 48; // chr(48) == '0'
	    pos++;
	}
	seqPos += seqCount;
    // misses
    while(MD[pos] == 'A' || MD[pos] == 'T' || MD[pos] == 'C' || MD[pos] == 'G' || MD[pos] == 'N'){
      double logMatch = loglikeMatch(readQual, seqPos, 1, qOff);
      double logMiss = 0.0;
      if(MD[pos] == 'A' || MD[pos] == 'T' || MD[pos] == 'C' || MD[pos] == 'G'){
        // correct likelihood for match in CIGAR
    	  logMiss = loglikeMiss(readQual, seqPos, 1, qOff);
      }
      else if(MD[pos] == 'N'){
    	  logMiss += log(0.25);
      } else {
    	  logMiss += log(0.25);
      }
      loglikelihood += logMiss - logMatch;
      seqPos++;
      pos++;
      //printf("MD %d miss  %d. %f\n", seqPos, 1, loglikelihood);
    }

    // deletions
    if(MD[pos] == '^'){
      pos++;
      while(isalpha(MD[pos])){
        pos++;
      }
    }

    // sees if we are at the end
    if(MD[pos] == '\0'){
      stop = 1;
    }
  }

  //printf("getMDLogLikelihood(%s): %e\n", MD, loglikelihood);
  // no assertion that logLikelihood is <= 0, as this is a correction to the logMatch already applied
  return loglikelihood;
}

double getCIGARLogLikelihoodBAM(int numCigarOperations, uint32_t *cigar, char *readQual, int qOff, int *inserts, int *deletions, int *totalMatch) {
  int i;
  int seqPos = 0;
  double logLikelihood = 0.0;
  for(i=0 ; i < numCigarOperations ; i++) {
    uint32_t cigarInt = *(cigar+i);
    uint32_t cigarFlag = (cigarInt & BAM_CIGAR_MASK);
    uint32_t count = (cigarInt >> BAM_CIGAR_SHIFT);
    //printf("CIGAR: %u %u %u\n", cigarInt, cigarFlag, count);
    switch (cigarFlag) {
      case(BAM_CMATCH) :
        *totalMatch += count;
        //likelihood *= likeMatch(readQual, seqPos, count, qOff);
        logLikelihood += loglikeMatch(readQual, seqPos, count, qOff);
        seqPos += count;
        break;
      case(BAM_CINS)   :
        *inserts += count;
        //likelihood *= likeInsertion(readQual, seqPos, count, qOff);
        logLikelihood += loglikeInsertion(readQual, seqPos, count, qOff);
        seqPos += count;
        break;
      case(BAM_CDEL)   :
        *deletions += count;
        //likelihood *= likeDeletion(readQual, seqPos, count, qOff);
        logLikelihood += loglikeDeletion(readQual, seqPos, count, qOff);
        // deletions do not increase seqPos
        break;
      case(BAM_CREF_SKIP):
        // assume this is a spliced alignment for RNA, so okay
        break;
      case(BAM_CHARD_CLIP):
        // clipped region is not in seq
        break;
      case(BAM_CSOFT_CLIP):
        //likelihood *= likeMiss(readQual, seqPos, count, qOff);
        logLikelihood += loglikeMiss(readQual, seqPos, count, qOff);
        seqPos += count;
        break;
    }
  }
  //double likelihood = exp(logLikelihood);
  //printf("getCIGARLikelihoodBAM(): %e, %e\n", likelihood, logLikelihood);
  assert(logLikelihood <= 0.0);
  return logLikelihood;
}

// takes in a read and returns the match loglikelihood (due to matches, mismatches, indels)
double getMatchLogLikelihoodBAM(bam1_t *read, int qOff){
  assert(read != NULL);
  double loglikelihood;
  // read CIGAR first
  char *readQual = (char*) bam1_qual(read);
  uint32_t *cigar = bam1_cigar(read);
  int inserts = 0;
  int deletions = 0;
  int totalMatch = 0;
  //printf("getMatchLikelihoodBAM(%s, %d)\n", bam1_qname(read), qOff);

  loglikelihood = getCIGARLogLikelihoodBAM(read->core.n_cigar, cigar, readQual, qOff, &inserts, &deletions, &totalMatch);
  assert(loglikelihood <= 0.0);

  char *md = (char*) bam_aux_get(read, "MD");
  //printf("%s %f MD:%s\n", bam1_qname(read), likelihood, md);
  if (md != NULL && md[0] == 'Z') {
	  loglikelihood += getMDLogLikelihood(md + 1, readQual, qOff);
  } else {
    printf("WARNING: could not find the MD tag for %s\n", bam1_qname(read));
  }

  //printf("getMatchLogLikelihoodBAM(%s, %d) = %e\n", bam1_qname(read), qOff, loglikelihood);
  return loglikelihood;
}

// returns the 2-bit hash representation of a nucl. given its place in the kmer
int kmerHash(char c1, int place){
  int hash = 0;
  switch(c1) {
  case 'A': case 'a': break;
  case 'T': case 't': hash = 0x1 << (place*2); break;
  case 'C': case 'c': hash = 0x2 << (place*2); break;
  case 'G': case 'g': hash = 0x3 << (place*2); break;
  default: hash = -1073741824;
  }
  return hash;
}

// builds up a 2*kmerLen-bit hash (for the seq) given a seqenence and starting position
int getKmerHash(char *seq, int startPos, int kmerLen){
  int i;
  int hash = 0;
  for(i = 0; i < kmerLen; i++){
    hash += kmerHash(seq[startPos+i], i);
  }
  return hash;
}

// calculates the kmer statistics
void computeKmerStats(assemblyT *theAssembly, int kmer){
  int i, j, k, totalKmers, hash;
  // calculate total possible kmers
  totalKmers = 1;
  for(i = 0; i < kmer; i++){
    totalKmers = totalKmers*4;
  }
  int *kmerVec = malloc(totalKmers*sizeof(int));
  // find all kmers present
  for(i = 0; i < theAssembly->numContigs; i++){
    contig_t *contig = theAssembly->contigs[i];
    // initialize kmerVec
    for(j = 0; j < totalKmers; j++){
      kmerVec[j] = 0;
    }
    totalKmers = 0;
    // add up the kmers
    for(j = 0; j < contig->seqLen - kmer + 1; j++){
      hash = getKmerHash(contig->seq, j, kmer);
      //printf("Hash = %i\n", hash);
      if(hash > -1){
        kmerVec[hash]++;
        totalKmers++;
      }
    }
    //printf("Calculated all %i kmers!\n", totalKmers);
    // calculate probability of seeing that kmer based on the rest of the contig
    // first kmer - 1 unrolled
    for(j = 0; j < kmer; j++){
      contig->kmerLikelihood[j] = 0.0;
      for(k = 0; k < j+1; k++){
        hash = getKmerHash(contig->seq, k, kmer);
        if(hash > -1){
          contig->kmerLikelihood[j] = contig->kmerLikelihood[j] + 1.0/(double)(j+1)*(double)(kmerVec[hash])/(double)(totalKmers);
        }
      }
      //printf("New likelihood[%i]: %f.\n", j, contig->kmerLikelihood[j]);
    }
    //printf("First.\n");
    // middle bunch
    for(j = kmer; j < contig->seqLen - kmer; j++){
      for(k = 0; k < kmer; k++){
        hash = getKmerHash(contig->seq, j - k, kmer);
        if(hash > -1){
          contig->kmerLikelihood[j] = contig->kmerLikelihood[j] + 1.0/(double)(kmer)*(double)(kmerVec[hash])/(double)(totalKmers);
        }
      }
      //printf("New likelihood[%i]: %f.\n", j, contig->kmerLikelihood[j]);
    }
    //printf("Mid.\n");
    // last bits
    for(j = contig->seqLen - kmer; j < contig->seqLen; j++){
      contig->kmerLikelihood[j] = 0.0;
      for(k = j - kmer + 1; k < j - kmer + 1 + (contig->seqLen - j); k++){
        hash = getKmerHash(contig->seq, k, kmer);
        if(hash > -1){
          contig->kmerLikelihood[j] = contig->kmerLikelihood[j] + 1.0/(double)(contig->seqLen - j)*(double)(kmerVec[hash])/(double)(totalKmers);
        }
      }
      //printf("New likelihood[%i]: %f.\n", j, contig->kmerLikelihood[j]);
    }
    //printf("Last.\n");
    totalKmers = 0;
  }
  free(kmerVec);
}

// used in getPlacementWinner
unsigned int JSHash(char* str)
{
  unsigned int hash = 1315423911;
  char c;

  while (str != NULL && (c=*str++) != '\0')
  {
    hash ^= ((hash << 5) + c + (hash >> 2));
  }

  return hash;
}

// finds the sum of the total likelihoods given the head of the list
double getTotalLikelihood(alignSet_t *head) {
  double likeNormalizer = 0.0;
  likeNormalizer += head->likelihood;
  alignSet_t *current = head;
  while(current->nextAlignment != NULL){
    current = current->nextAlignment;
    likeNormalizer += current->likelihood;
  }
  //printf("Normalizer: %f\n", likeNormalizer);
  return likeNormalizer;
}

// find the placement winner
alignSet_t *getPlacementWinner(alignSet_t *head, double likeNormalizer, int *winner) {
  assert(head != NULL && winner != NULL);

  *winner = -1;
  if(likeNormalizer == 0.0){ // no real placement
    return NULL;
  }

  // instead of purely random choice, choose a consistent, but still random choice
  unsigned int iseed = JSHash(head->name);
  srand(iseed); // instead of: srand ((unsigned int) time(NULL));

  double tRand = ( likeNormalizer * rand() / ( RAND_MAX + 1.0 ) );
  double soFar = 0.0;

  int i = 0;
  alignSet_t *current = head;
  if(head->likelihood > tRand){
    *winner = 0;
  }else{
    soFar += head->likelihood;
    while(current->nextAlignment != NULL){
      current = current->nextAlignment;
      i++;
      if(current->likelihood + soFar > tRand){
        *winner = i;
        break;
      }else{
        soFar += current->likelihood;
      }
    }
  }
  if(current->likelihood >= 0.0){
    return current;
  }else{ // not a real placement
    return NULL;
  }
}

// apply statistics
void applyDepthAndMatchToContig(alignSet_t *alignment, assemblyT *theAssembly, double likeNormalizer) {
  double likelihood = alignment->likelihood;
  assert(likelihood >= 0);
  int j;
  if (alignment->start1 >= 0) {
    assert(alignment->contigId1 >= 0 && alignment->contigId1 < theAssembly->numContigs);
    contig_t *contig1 = theAssembly->contigs[alignment->contigId1];
    assert(alignment->start1 < contig1->seqLen);
    assert(alignment->end1 <= contig1->seqLen);
    assert(alignment->start1 < alignment->end1);
    for(j = alignment->start1; j < alignment->end1; j++){
      contig1->depth[j] += 1.0; // We picked a winner, it gets full prob
      contig1->matchLikelihood[j] += log(likelihood);
    }
  }
  if (alignment->start2 >= 0) {
	assert(alignment->contigId2 >= 0 && alignment->contigId2 < theAssembly->numContigs);
	contig_t *contig2 = theAssembly->contigs[alignment->contigId2];
	assert(alignment->start1 < contig2->seqLen);
	assert(alignment->end1 <= contig2->seqLen);
	assert(alignment->start1 < alignment->end1);
    for(j = alignment->start2; j < alignment->end2; j++){
      contig2->depth[j] += 1.0;
      contig2->matchLikelihood[j] += log(likelihood);
    }
  }
}

// this applies the placement(s) to the assembly part (SINGLE PART)
// I feel like this could be sped up with a hash table vs the current linked lists, but we will see...
int applyPlacement(alignSet_t *head, assemblyT *theAssembly){

  // normalize the probs
  double likeNormalizer = getTotalLikelihood(head);

  int winner = -1;
  alignSet_t *current = getPlacementWinner(head, likeNormalizer, &winner);

  if(current == NULL){
    printf("No winner, failed to place %s. currentLikelihood: %f, Normalizer: %f\n", head->name, head->likelihood, likeNormalizer);
    return -1;
  }

  // apply the placement
  applyDepthAndMatchToContig(current, theAssembly, likeNormalizer);

  return winner;
}

// compute the depth statistics
int computeDepthStats(assemblyT *theAssembly){
  int i, j;
  double depthNormalizer[102];
  long depthNormalizerCount[102];
  for(i = 0; i < 101; i++){
    depthNormalizer[i] = 0.0;
    depthNormalizerCount[i] = 0;
  }
  double tempLike;
  long tooLowCoverageBases = 0;
  long noGCInformation = 0;
  for(i = 0; i < theAssembly->numContigs; i++){ // for each contig
    contig_t *contig = theAssembly->contigs[i];
    for(j = 0; j < contig->seqLen; j++){
      float depth = contig->depth[j];
      if (depth < 0.1) {
    	tooLowCoverageBases++;
        continue;
      }
      int GCpct = contig->GCcont[j];
      if (GCpct > 100) {
    	  noGCInformation++;
    	  continue;
      }
      depthNormalizer[GCpct] += depth;
      depthNormalizerCount[GCpct] += 1;
    }
  }
  for(j = 0; j < 101; j++){
    if(depthNormalizerCount[j] > 0){
      depthNormalizer[j] = depthNormalizer[j]/(double)depthNormalizerCount[j];
    }else{
      depthNormalizer[j] = 0.1;
    }
    printf("depth at GC[%d] = %f (%ld samples)\n", j, depthNormalizer[j], depthNormalizerCount[j]);
  }
  for(i = 0; i < theAssembly->numContigs; i++){ // for each contig
    contig_t *contig = theAssembly->contigs[i];
    for(j = 0; j < contig->seqLen; j++){
      int GCpct = contig->GCcont[j];
      if (GCpct > 100)
    	  continue;
      tempLike = poissonPMF(contig->depth[j], depthNormalizer[GCpct]); // log poisson pmf
      if(tempLike < minLogLike || isnan(tempLike)){tempLike = minLogLike;}
      contig->depthLikelihood[j] = tempLike;
      tempLike = contig->matchLikelihood[j]/contig->depth[j]; // log applied in applyDepthAndMatchToContig()
      if(tempLike < minLogLike || isnan(tempLike)){tempLike = minLogLike;}
      contig->matchLikelihood[j] = tempLike;
    }
  }
  printf("bases with too low coverage: %ld\n", tooLowCoverageBases);
  printf("bases with no GC metric (small contigs): %ld\n", noGCInformation);
  return 1;
}

int guessQualityOffset(bam1_t *read) {
  int len = read->core.l_qseq;
  if (len <= 0)
    return -1;
  int qualOffset = -1;
  const int maxExpectedQuality = 50;
  char *qualSeq = (char*) bam1_qual(read);
  int i;
  for(i=0; i < len; i++) {
    if (qualSeq[i] > 64 && qualSeq[i] > 33+maxExpectedQuality)
      qualOffset = 64;
    else if (qualSeq[i] > 0 && qualSeq[i] < 33)
      qualOffset = 0;
    else if (qualSeq[i] > 33 && qualSeq[i] < maxExpectedQuality+33 && qualSeq[i] < 64 && qualSeq[i] > maxExpectedQuality+0)
      qualOffset = 33;

    if (qualOffset >= 0) {
      printf("guessed quality offset is %d\n", qualOffset);
      break;
    }
  }
  return qualOffset;
}

libraryParametersT *computeLibraryParameters(samfile_t *ins, double outlierFraction, int qOff) {

  int i,j;
  long *mapLens[MATE_ORIENTATION_MAX];

  // allocate memory
  for(i=0; i < MATE_ORIENTATION_MAX; i++)
    mapLens[i] = malloc(sizeof(long)*mapLens_MAX);
  libraryParametersT *libParams = malloc(sizeof(libraryParametersT));

  // initialize
  for(j=0; j < MATE_ORIENTATION_MAX; j++) {
    for(i = 0; i < mapLens_MAX; i++)
      mapLens[j][i] = 0;
    libParams->mateParameters[j].insertLength = 0.0;
    libParams->mateParameters[j].insertStd = 0.0;
    libParams->mateParameters[j].count = 0;
    libParams->mateParameters[j].isValid = 0;
  }
  libParams->qOff = qOff;
  libParams->avgReadSize = 0;
  libParams->numReads = 0;
  libParams->isSortedByName = -1; // undefined

  bam1_t *thisRead = bam_init1();
  bam1_t *thisReadMate = bam_init1();

  long readCount = 0;
  long improperPair = 0;

  while(1){
    enum MATE_ORIENTATION orientation = readMatesBAM(ins, libParams, thisRead, thisReadMate);
    if (orientation == NO_READS)
      break;
    else if (orientation == UNMAPPED_PAIR)
      continue;

    libraryMateParametersT *mateParams = &libParams->mateParameters[orientation];

    int mapLen = -1;
    int numReads = 0;
    switch(orientation) {
      case(VALID_FR):
      case(VALID_RF):
      case(VALID_FF):
      case(NOT_PROPER_FR):
      case(NOT_PROPER_RF):
      case(NOT_PROPER_FF):
        if (libParams->isSortedByName == -1) {
          libParams->isSortedByName = 1;
          printf("Setting library to be sorted by name\n");
        }
        mapLen = getMapLenBAM(thisRead);
        numReads = 2;
        break;
      case(HALF_VALID_MATE):
        numReads = 1;
        if (libParams->isSortedByName == -1) {
          libParams->isSortedByName = 0;
          printf("Setting library to be unsorted by name\n");
        }
        break;
      case(UNRELATED_PAIR):
        numReads = 2;
        if (libParams->isSortedByName == -1) {
          libParams->isSortedByName = 0;
          printf("Setting library to be unsorted by name\n");
        }
        break;
      case(READ1_ONLY):
      case(READ2_ONLY):
      case(SINGLE_READ):
      case(UNMAPPED_SINGLE):
        numReads = 1;
        break;
      case(CHIMER):
      case(UNMAPPED_PAIR):
      default:
    	numReads = 2;
        improperPair++;
    }

    libParams->avgReadSize += getSeqLenBAM(thisRead) + getSeqLenBAM(thisReadMate);
    libParams->numReads += numReads;
    if (libParams->qOff < 0) {
      libParams->qOff = guessQualityOffset(thisRead);
    }
    if (libParams->qOff < 0) {
      libParams->qOff = guessQualityOffset(thisReadMate);
    }

    mateParams->count++;
    if (mapLen > 0) {
      mateParams->insertLength += mapLen;
      if (mapLen < mapLens_MAX)
        ++mapLens[orientation][mapLen];
    }

    if ((++readCount & 0xfffff) == 0)
      printf("Read %ld reads...\n", readCount);
  }
  bam_destroy1(thisRead);
  bam_destroy1(thisReadMate);

  // zero out top and bottom outliers
  double totalValid = 0.0;
  for(j = 0; j < MATE_ORIENTATION_MAX; j++) {
    libraryMateParametersT *mateParams = &libParams->mateParameters[j];
    printf("Evaluating %s orientation with %ld reads\n", MATE_ORIENTATION_LABELS[j], mateParams->count);
    long observed = 0;
    long purged = 0;
    long lengthTotal = 0;
    for(i = 0; i < mapLens_MAX; i++) {
      if (observed > (1.0-outlierFraction) * mateParams->count) {
        purged += mapLens[j][i];
        mapLens[j][i] = 0;
      } else {
        observed += mapLens[j][i];
        lengthTotal += i*mapLens[j][i];
      }
      if (observed < outlierFraction * mateParams->count) {
        purged += mapLens[j][i];
        mapLens[j][i] = 0;
      }
    }
    mateParams->libraryFraction = (double) mateParams->count * 2.0 / (double) libParams->numReads;
    if (j == READ1_ONLY || j == READ2_ONLY || j == SINGLE_READ)
      mateParams->libraryFraction /= 2.0;

    if (mateParams->count > 0 && ( j == SINGLE_READ || j == VALID_FR || j == VALID_RF || j == VALID_FF || j == NOT_PROPER_FR || j == NOT_PROPER_RF || j == NOT_PROPER_FF )) {

      // TODO better test for significance and normal distribution for a valid orientation
      if (mateParams->libraryFraction > SIGNIFICANT_LIBRARY_FRACTION) {
        totalValid += mateParams->count * (j == SINGLE_READ ? 1.0 : 2.0);
        mateParams->isValid = 1;
      } else {
    	mateParams->isValid = 0;
      }

      printf("Read %ld properly oriented sequences (%s) and purged %ld %0.1lf%% & %0.1lf%% outliers.  This %s a valid orientation (%01lf%%)\n",
          mateParams->count,
          MATE_ORIENTATION_LABELS[j],
          purged,
          (outlierFraction*100.0),
          ((1.0-outlierFraction)*100.0),
          mateParams->isValid ? "is" : "is NOT",
          mateParams->libraryFraction*100.0);
    }

    long modifiedReadCount = mateParams->count - purged;

    if (mateParams->isValid) {
      mateParams->insertLength = (double)lengthTotal/(double)modifiedReadCount;
      for(i = 0; i < mapLens_MAX; i++){
        if(mapLens[j][i] > 0){
          double tmp = mapLens[j][i]*((double)i - mateParams->insertLength)*((double)i - mateParams->insertLength);
          mateParams->insertStd += tmp;
          printf("i : %s mapLens[i] :: %i : %ld\n", MATE_ORIENTATION_LABELS[j], i, mapLens[j][i]);
        }
      }
      mateParams->insertStd = sqrt(mateParams->insertStd/(double)(modifiedReadCount-1));
      printf("Found %s sample avg insert length to be %lf from %ld mapped reads\n", MATE_ORIENTATION_LABELS[j], mateParams->insertLength, modifiedReadCount);
      printf("Found %s sample insert length std to be %lf\n", MATE_ORIENTATION_LABELS[j], mateParams->insertStd);
    }
  }
  printf("There were %ld improper pairs and unmapped reads\n", improperPair);
  libParams->avgReadSize = libParams->avgReadSize / libParams->numReads;
  printf("Found sample avg read size to be %ld\n", libParams->avgReadSize);

  double totalSingle = libParams->mateParameters[READ1_ONLY].count + libParams->mateParameters[READ2_ONLY].count;
  double totalUnmapped = totalSingle;
  totalSingle += libParams->mateParameters[SINGLE_READ].count;
  totalUnmapped += libParams->mateParameters[UNMAPPED_SINGLE].count + libParams->mateParameters[UNMAPPED_PAIR].count;

  libParams->totalValidSingleFraction = totalSingle / libParams->numReads;
  libParams->totalValidMateFraction = (totalValid - totalSingle) / libParams->numReads;
  libParams->totalUnmappedFraction = totalUnmapped / libParams->numReads;
  libParams->totalChimerMateFraction = 1.0 - libParams->totalValidMateFraction - libParams->totalValidSingleFraction - libParams->totalUnmappedFraction;

  printf("ValidMates %0.3lf%%, Single %0.3lf%%, ChimerMates %0.3lf%%, Unmapped %0.3lf%%\n",
		  libParams->totalValidMateFraction*100,
		  libParams->totalValidSingleFraction*100,
		  libParams->totalChimerMateFraction*100,
		  libParams->totalUnmappedFraction*100);

  // release memory
  for(i=0; i < MATE_ORIENTATION_MAX; i++)
    free(mapLens[i]);

  return libParams;
}

int mateTreeCmp(const void *pa, const void *pb) {
  assert(pa != NULL && pb != NULL);
  alignSet_t *a = (alignSet_t*) pa;
  alignSet_t *b = (alignSet_t*) pb;
  return strcmp(a->name, b->name);
}

  alignSet_t *getOrStoreMateAlignment(void **mateTree, alignSet_t *thisAlignment, bam1_t *thisRead) {
    if (mateTree == NULL)
      return NULL;
    assert(thisAlignment->start2 == -1);
    void *found = NULL;
    alignSet_t *stored = NULL;
    if (thisRead->core.pos >= thisRead->core.mpos) {
      found = (alignSet_t*) tfind((void*)thisAlignment, mateTree, mateTreeCmp);
    }
    if (found == NULL) {
      stored = (alignSet_t*) malloc(sizeof(alignSet_t));
      if (stored == NULL) {
        printf("Unable to malloc\n");
        exit(1);
      }
      copyAlignment(stored, thisAlignment);
      void *successful = tsearch((void*)stored, mateTree, mateTreeCmp);
      assert(successful != NULL);
      mateTreeCount++;
      //printf("Stored %s\n", stored->name);
      return NULL; // i.e. HALF_VALID_MATE;
    } else {
      alignSet_t *mateAlignment = *((alignSet_t**) found);
      //printf("Found %s %s\n", thisAlignment->name, mateAlignment->name);
      assert(strcmp(mateAlignment->name, thisAlignment->name) == 0);
      // get the read2 parameters (from read1 in found)
      thisAlignment->start2 = mateAlignment->start1;
      thisAlignment->end2 = mateAlignment->end2;
      tdelete((void*)mateAlignment, mateTree, mateTreeCmp);
      return mateAlignment;
      mateTreeCount--;
    }
  }

void mateTreeFreeNode(void *nodep) {
  if (nodep != NULL) {
    mateTreeCount--;
    //printf("mateTreeFreeNode(%p) freeing (%d)\n", nodep, mateTreeCount);
    free(nodep);
  }
}

void setSingleRead2Alignment(bam_header_t *header, alignSet_t *read2Only, alignSet_t *thisAlignment, bam1_t *thisReadMate, double likelihood) {
  // transfer thisAlignment read2 to read2Only read1
  assert(thisReadMate->core.tid == thisAlignment->contigId2);
  read2Only->start1 = thisAlignment->start2;
  read2Only->end1 = thisAlignment->end2;
  read2Only->contigId1 = thisAlignment->contigId2;

  // reset both alignments read2 coordinates
  thisAlignment->start2 = thisAlignment->end2 = thisAlignment->contigId2 = -1;
  read2Only->start2 = read2Only->end2 = read2Only->contigId2 = -1;

  read2Only->likelihood = likelihood;
  strncpy(read2Only->name, bam1_qname(thisReadMate), MAX_NAME_LENGTH);
  read2Only->name[MAX_NAME_LENGTH] = '\0'; // ensure null termination

  read2Only->nextAlignment = NULL;
}

enum MATE_ORIENTATION setAlignment(bam_header_t *header, assemblyT *theAssembly, alignSet_t *thisAlignment, alignSet_t *secondaryAlignment, void **mateTree, libraryParametersT *libParams, enum MATE_ORIENTATION orientation, bam1_t *thisRead, bam1_t *thisReadMate) {
  assert(thisAlignment != NULL);
  double loglikelihoodRead1 = 0.0;
  double loglikelihoodRead2 = 0.0;

  double likelihoodInsert;

  strncpy(thisAlignment->name, bam1_qname(thisRead), MAX_NAME_LENGTH);
  thisAlignment->name[MAX_NAME_LENGTH] = '\0'; // ensure null termination
  thisAlignment->nextAlignment = NULL;

  // reset any existing secondaryAlignment
  secondaryAlignment->likelihood = 0.0;

  int qOff = libParams->qOff;
  thisAlignment->likelihood = 1.0;
  if (thisRead == NULL || (thisRead->core.flag & BAM_FUNMAP) == BAM_FUNMAP) {
    thisAlignment->start1 = -1;
    thisAlignment->end1   = -1;
  } else {
    loglikelihoodRead1  = getMatchLogLikelihoodBAM(thisRead, qOff);
    thisAlignment->start1 = thisRead->core.pos;
    thisAlignment->end1   = bam_calend(&thisRead->core, bam1_cigar(thisRead));
    assert(thisAlignment->start1 <= thisAlignment->end1);
    thisAlignment->contigId1 = thisRead->core.tid;
  }

  if (thisReadMate == NULL || (thisReadMate->core.flag & BAM_FUNMAP) == BAM_FUNMAP) {
    thisAlignment->start2 = -1;
    thisAlignment->end2   = -1;
  } else {
    loglikelihoodRead2  = getMatchLogLikelihoodBAM(thisReadMate, qOff);
    thisAlignment->start2 = thisReadMate->core.pos;
    thisAlignment->end2   =  bam_calend(&thisReadMate->core, bam1_cigar(thisReadMate));
    assert(thisAlignment->start2 <= thisAlignment->end2);
    if (thisAlignment->start1 < 0) {
      thisAlignment->contigId2 = thisReadMate->core.tid;
    }
  }

  libraryMateParametersT *mateParameters = &libParams->mateParameters[orientation];
  double logzNormalizeRead1, logzNormalizeRead2, logzNormalizeReadMates;

  switch (orientation) {
    case (VALID_FR):
    case (VALID_RF):
    case (VALID_FF):
    case (NOT_PROPER_FR):
    case (NOT_PROPER_RF):
    case (NOT_PROPER_FF):
      // two reads

      if (mateParameters->isValid) {
        // valid orientation
        logzNormalizeReadMates = logzNormalizationReadQual(thisRead, thisReadMate, libParams->qOff);
        likelihoodInsert = getInsertLikelihoodBAM(thisRead, mateParameters->insertLength, mateParameters->insertStd);
        thisAlignment->likelihood = libParams->totalValidMateFraction * likelihoodInsert;
        thisAlignment->likelihood *= exp( loglikelihoodRead1 + loglikelihoodRead2 - logzNormalizeReadMates);
        break;
      } else {
        // change the orientation... this is actually a chimer
        orientation = CHIMER;
      }
      //  continue... this is actually a chimer
    case (CHIMER) :
      //printf("WARNING: chimeric read mate pair %s.\n", bam1_qname(thisRead));

      likelihoodInsert = libParams->totalChimerMateFraction;
      if (thisRead->core.tid == thisRead->core.mtid && theAssembly->contigs[thisRead->core.tid]->isCircular == 1) {
        // TODO refine based on proximity to end of contigs in a circular genome
        // i.e. factor in: likelihoodThatRead1AndRead2AreCloseToContigEdgeSoAreNotChimersButProbablySpanningMatePairs
      }
      logzNormalizeRead1 = logzNormalizationReadQual(thisRead, NULL, libParams->qOff);
      thisAlignment->likelihood *= likelihoodInsert * exp(loglikelihoodRead1 - logzNormalizeRead1);

      // set secondaryAlignment to map both reads separately...
      if (thisReadMate != NULL) {
          logzNormalizeRead2 = logzNormalizationReadQual(thisReadMate, NULL, libParams->qOff);
          double likelihoodMate = likelihoodInsert * exp(loglikelihoodRead2 - logzNormalizeRead2);
          setSingleRead2Alignment(header, secondaryAlignment, thisAlignment, thisReadMate, likelihoodMate);
      }

      break;

    case (UNRELATED_PAIR):
    case (READ1_ONLY):
    case (READ2_ONLY):

      if(thisReadMate != NULL && (thisReadMate->core.flag & BAM_FUNMAP) != BAM_FUNMAP) {
        // set seconaryAlignment to thisReadMate, remove thisReadMate from thisAlignment
    	double logzNormalizeRead2 = logzNormalizationReadQual(thisReadMate, NULL, libParams->qOff);
        setSingleRead2Alignment(header, secondaryAlignment, thisAlignment, thisReadMate, exp(loglikelihoodRead2 - logzNormalizeRead2));

        if ((thisReadMate->core.flag & BAM_FMUNMAP) == BAM_FMUNMAP) {
          // mate is not mapped, apply likelihood now
          secondaryAlignment->likelihood *= libParams->totalValidSingleFraction;
        } else {
          //assert(thisRead == NULL || thisReadMate->core.tid == thisRead->core.mtid);
          // store this or get mate if already seen
          alignSet_t *mateAlignment = getOrStoreMateAlignment(mateTree, secondaryAlignment, thisReadMate);
          if (mateAlignment != NULL) {
            double matelikelihood = mateAlignment->likelihood;
            free(mateAlignment);
            enum MATE_ORIENTATION thisOrientation = getPairedMateOrientation(thisReadMate);
            mateParameters = &libParams->mateParameters[thisOrientation];
            likelihoodInsert = getInsertLikelihoodBAM(thisReadMate, mateParameters->insertLength, mateParameters->insertStd);
            secondaryAlignment->likelihood *= matelikelihood * likelihoodInsert;
          } else {
            // do not process this one yet
            secondaryAlignment->likelihood = 0.0;
          }
        }
      }

      if (thisRead != NULL && (thisRead->core.flag & BAM_FUNMAP) != BAM_FUNMAP) {
    	logzNormalizeRead1 = logzNormalizationReadQual(thisRead, NULL, libParams->qOff);
        thisAlignment->likelihood = exp(loglikelihoodRead1 - logzNormalizeRead1);
        thisAlignment->start2 = thisAlignment->end2 = -1;
        if ((thisRead->core.flag & BAM_FMUNMAP) == BAM_FMUNMAP) {
          // mate is not mapped, apply likelihood now
          likelihoodInsert = libParams->totalValidSingleFraction;
          thisAlignment->likelihood *= likelihoodInsert;
        } else {
          //assert(thisReadMate == NULL || thisRead->core.tid == thisReadMate->core.mtid);
          // store this or get mate if already seen
          alignSet_t *mateAlignment = getOrStoreMateAlignment(mateTree, thisAlignment, thisRead);
          if (mateAlignment != NULL) {
            double likelihoodRead2 = mateAlignment->likelihood;
            free(mateAlignment);

            // reset orientation, if needed
            orientation = getPairedMateOrientation(thisRead);
            mateParameters = &libParams->mateParameters[orientation];
            likelihoodInsert = getInsertLikelihoodBAM(thisRead, mateParameters->insertLength, mateParameters->insertStd);
            thisAlignment->likelihood *= likelihoodRead2 * likelihoodInsert;
          } else {
            // do not process this one yet
            thisAlignment->likelihood = 0.0;
            orientation = HALF_VALID_MATE;
          }
        }
      }
      break;

    case (SINGLE_READ):
      logzNormalizeRead1 = logzNormalizationReadQual(thisRead, NULL, libParams->qOff);
      thisAlignment->likelihood = libParams->totalValidSingleFraction * exp(loglikelihoodRead1 - logzNormalizeRead1);
      if (thisReadMate != NULL) {
    	logzNormalizeRead2 = logzNormalizationReadQual(thisReadMate, NULL, libParams->qOff);
        setSingleRead2Alignment(header, secondaryAlignment, thisAlignment, thisReadMate, libParams->totalValidSingleFraction * exp(loglikelihoodRead2 - logzNormalizeRead2));
      }
      break;

    case (UNMAPPED_PAIR):
      break;

    default :
      thisAlignment->likelihood = 0.0;
      if (thisRead == NULL) { printf("thisread is null!!! Skipping %s\n", MATE_ORIENTATION_LABELS[orientation]); break; }
      if (thisReadMate == NULL) { break; }
      assert(thisRead != NULL && thisReadMate != NULL);
      printf("Skipping %s read %s %s\n", MATE_ORIENTATION_LABELS[orientation], bam1_qname(thisRead), bam1_qname(thisReadMate));
      break;
  }

  return orientation;
}

// NORMALIZATION
// divide by the expected loglikelihood of the read by the normalization factor Z (from Bayes rule)
// given only its length and the parameters of the distributions (See paper appendix)
double logzNormalizationReadQual(bam1_t *thisRead, bam1_t *thisReadMate, int qOff){
  // find the average quality to save computation/precision in combinatorics
  double Qavg = 0.0;
  char *readQual = (char*) bam1_qual(thisRead);
  int totalLen = thisRead->core.l_qseq;
  int i;
  for(i = 0; i < thisRead->core.l_qseq; i++){
    Qavg += getQtoP(readQual[i], qOff);
  }
  if(thisReadMate != NULL){
    readQual = (char*) bam1_qual(thisReadMate);
    totalLen += thisReadMate->core.l_qseq;
    for(i = 0; i < thisReadMate->core.l_qseq; i++){
      Qavg += getQtoP(readQual[i], qOff);
    }
  }
  Qavg = Qavg/(double)totalLen;
  double QmisMatch = (1.0 - Qavg)/3.0;// assume all remaining probability goes equally to the remaining bases

  // find the expected match score
  //double expMatch = 0.0;
  //for(i = 0; i < totalLen; i++){
  //  expMatch += pow(pow(Qavg, i)*pow(QmisMatch, totalLen - i - 1), 2);
  //}

  //printf("Qavg: %lf %lf\n", Qavg, QmisMatch);
  //printf("expMatch: %e %e\n", expMatch, maxExpMatch);

  double logQavg = log(Qavg);
  double logQmisMatch = log(QmisMatch);

  // normalize over the maximum value to prevent double precision underflows
  //double maxExpMatch = pow(pow(Qavg > QmisMatch ? Qavg : QmisMatch, totalLen - 1)*1,2);
  double logMaxExpMatch = 2.0*(totalLen-1)*(logQavg > logQmisMatch ? logQavg : logQmisMatch);

  // find the log expected match score
  double tmpExpMatch = 0.0;
  for(i = 0; i < totalLen; i++){
	  //log(pow(pow(Qavg, i)*pow(QmisMatch, totalLen - i - 1), 2));
	  double logTmp = 2.0*(i*logQavg + (totalLen-i-1)*logQmisMatch);
	  tmpExpMatch += exp( logTmp - logMaxExpMatch );
  }
  double logExpMatch = logMaxExpMatch + log(tmpExpMatch);

  //printf("logQavg: %lf %lf\n", logQavg, logQmisMatch);
  //printf("logExpMatch: %lf %lf\n",  logExpMatch, logMaxExpMatch);
  //printf("expMatch: %e %e\n", exp(logExpMatch), exp(logMaxExpMatch));
  //if (exp(logExpMatch) > expMatch * 1.00001 || exp(logExpMatch) < expMatch * 0.999999)
  //	  printf("expMatch: %e %e %e\n", expMatch - exp(logExpMatch), expMatch, exp(logExpMatch));
  return logExpMatch;
}

double zNormalizationInsertStd(libraryMateParametersT *mateParams) {
  // int((pmf of normal(0,sigma))^2, 0, infty)
  double expIns = 1.0;
  if (mateParams->insertStd > 0)
	  expIns /= (2*sqrt(3.14159265)*mateParams->insertStd);

  return expIns;
}

void computeReadPlacements(samfile_t *ins, assemblyT *theAssembly, libraryParametersT *libParams, samfile_t *placementBam) {
  // initialize variables
  int i;
  alignSet_t alignments[N_PLACEMENTS];
  bam1_t *samReadPairs[N_PLACEMENTS*2];
  for(i=0; i < N_PLACEMENTS; i++) {
    samReadPairs[i*2] = bam_init1();
    samReadPairs[i*2+1] = bam_init1();
    initAlignment(&alignments[i]);
  }
  int samReadPairIdx = 0;

  alignSet_t *currentAlignment = NULL;
  alignSet_t *head = currentAlignment;
  alignSet_t secondaryAlignment;
  initAlignment(&secondaryAlignment);

  int failedToPlace = 0;
  int placed = 0;

  void *mateTree = NULL;

  int readCount = 0;
  while(1){
    bam1_t *thisRead = samReadPairs[samReadPairIdx*2];
    bam1_t *thisReadMate = samReadPairs[samReadPairIdx*2+1];
    alignSet_t *thisAlignment = &alignments[samReadPairIdx];
    samReadPairIdx++;

    enum MATE_ORIENTATION orientation = readMatesBAM(ins, libParams, thisRead, thisReadMate);
    if ((++readCount & 0xfffff) == 0)
      printf("Read %d reads...\n", readCount);
    if (orientation == NO_READS)
      break;

    orientation = setAlignment(ins->header, theAssembly, thisAlignment, &secondaryAlignment, &mateTree, libParams, orientation, thisRead, libParams->isSortedByName == 1 ? thisReadMate : NULL);
    if (orientation == UNMAPPED_PAIR) {
      samReadPairIdx--;
      continue;
    }

    // ************* //
    // NORMALIZATION //
    // ************* //

    // divide by the expected likelihood of the read by the normalization factor Z (from Bayes rule)
    // (read-qual normalization happens in setAlignment(q)
    // given only its length and the parameters of the distributions (See paper appendix)
    thisAlignment->likelihood /= zNormalizationInsertStd(&libParams->mateParameters[orientation]);
    

    if (orientation == NO_READS)
      break;

    //printf("%s\n", bam1_qname(thisRead));

    //process secondaryAlignment (thisReadMate separately)
    // apply placement of read2 to target2...
    if (thisReadMate != NULL && secondaryAlignment.likelihood > 0.0) {
      int winner = applyPlacement(&secondaryAlignment, theAssembly);
      if (winner < 0) {
        printf("WARNING: no placement found for read2 of chimer %s!\n", secondaryAlignment.name);
      } else {
        if (placementBam != NULL){
          bam_write1(placementBam->x.bam, thisReadMate);
        }
      }
    }

    if (orientation == HALF_VALID_MATE || thisAlignment->likelihood == 0.0) {
      // do not bother placing, just read the next one.
      samReadPairIdx--;
      continue;
    }

    //printf("Likelihoods (%s): %12f %12f %12f\n", bam1_qname(thisRead), likelihoodRead1, likelihoodRead2, likelihoodInsert);
    //printf("%s : %s .\n", currentAlignment->name, read.readName);

    // organize linked list of alignments based on current state of input stream
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
          if (alignments[winner].start1 >= 0)
            bam_write1(placementBam->x.bam, samReadPairs[winner*2]);
          if (alignments[winner].start2 >= 0)
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
    }
  }

  printf("%i maps placed, %i maps failed to place.\n", placed, failedToPlace);
  // tear down SAM/BAM variables
  for(i=0; i < N_PLACEMENTS; i++) {
    if (samReadPairs[i*2] != NULL)
      bam_destroy1(samReadPairs[i*2]);
    if (samReadPairs[i*2+1] != NULL)
      bam_destroy1(samReadPairs[i*2+1]);
  }

  if (mateTreeCount > 0) {
    printf("Applying placement for remaining/missing mate pairs (%d)\n", mateTreeCount);
    // TODO twalk stragglers, apply their placement as singles (should not get here on proper BAMs)
  }
  tdestroy(mateTree,mateTreeFreeNode);
  printf("Destroyed mateTree (%d)\n", mateTreeCount);
  //assert(mateTreeCount == 0);
}
