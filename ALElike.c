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

// finds the poisson pmf at value k for mean lambda
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
  return likelihood;
}

// finds the likelihood of a string of misses in read from seqPos to matchLen
double likeMiss(char *readQual, int seqPos, int missLen, int qOff){
  int i;
  double likelihood = 1.0;
  for(i = seqPos; i < seqPos + missLen; i++){
    //likelihood = likelihood*((1.0 - 0.99)/3.0);
    likelihood = likelihood*((1.0 - getQtoP(readQual[i], qOff))/3.0); // sometimes want -33/64?
  }
  return likelihood;
}

// finds the likelihood of a string of matches in read from seqPos to matchLen
double likeMatch(char *readQual, int seqPos, int matchLen, int qOff){
  int i;
  double likelihood = 1.0;
  for(i = seqPos; i < seqPos + matchLen; i++){
    //likelihood = likelihood*(0.99);
    likelihood = likelihood*(getQtoP(readQual[i], qOff));// sometimes want -33/64?
    //printf("likeMatch %d %c %f %f\n", i, readQual[i], QtoP[readQual[i] - qOff], likelihood);
  }
  //printf("likeMatch: %s %lf\n", read->readName, likelihood);
  return likelihood;
}

// finds the likelihood of an insertion (right now it is the same as a miss)
double likeInsertion(char *readQual, int seqPos, int insertionLength, int qOff) {
  // assume as unlikely as a substitution
  // TODO refine
  return likeMiss(readQual, seqPos, insertionLength, qOff);
}

// finds the likelihood of an deletion (right now it is the same as a miss)
double likeDeletion(char *readQual, int seqPos, int deletionLength, int qOff) {
  // assume as unlikely as a substitution of previous base
  // TODO refine
  assert(seqPos > 0);
  return pow(likeMiss(readQual, seqPos - 1, 1, qOff), (double)deletionLength);
}

// used to reduce likelihood in case of missmatches only
// (CIGAR already has accounted for matchlength, inserts, deletions)
double getMDLikelihood(char *MD, char *readQual, int qOff) {
  assert(MD != NULL && MD[0] != '\0');
  assert(readQual != NULL && readQual[0] != '\0');

  int stop = 0;
  int pos = 0;
  int seqPos = 0;
  double likelihood = 1.0;

  // parse MD field
  while(stop == 0){
    // matches
    while(isdigit(MD[pos])){
      pos++;
      seqPos++;
    }
    // misses
    while(MD[pos] == 'A' || MD[pos] == 'T' || MD[pos] == 'C' || MD[pos] == 'G' || MD[pos] == 'N'){
      if(MD[pos] == 'A' || MD[pos] == 'T' || MD[pos] == 'C' || MD[pos] == 'G'){
        // correct likelihood for match in CIGAR
        likelihood = likelihood * likeMiss(readQual, seqPos, 1, qOff) / likeMatch(readQual, seqPos, 1, qOff);
      }
      if(MD[pos] == 'N'){
        likelihood = likelihood*0.25 / likeMatch(readQual, seqPos, 1, qOff);
      }
      seqPos++;
      pos++;
      //printf("MD %d miss  %d. %f\n", seqPos, 1, likelihood);
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

  return likelihood;
}

double getCIGARLikelihoodBAM(int numCigarOperations, uint32_t *cigar, char *readQual, int qOff, int *inserts, int *deletions, int *totalMatch) {
  int i;
  int seqPos = 0;
  double likelihood = 1.0;
  for(i=0 ; i < numCigarOperations ; i++) {
    uint32_t cigarInt = *(cigar+i);
    uint32_t cigarFlag = (cigarInt & BAM_CIGAR_MASK);
    uint32_t count = (cigarInt >> BAM_CIGAR_SHIFT);
    //printf("CIGAR: %u %u %u\n", cigarInt, cigarFlag, count);
    switch (cigarFlag) {
      case(BAM_CMATCH) :
        *totalMatch += count;
        likelihood *= likeMatch(readQual, seqPos, count, qOff);
        seqPos += count;
        break;
      case(BAM_CINS)   :
        *inserts += count;
        likelihood *= likeInsertion(readQual, seqPos, count, qOff);
        seqPos += count;
        break;
      case(BAM_CDEL)   :
        *deletions += count;
        likelihood *= likeDeletion(readQual, seqPos, count, qOff);
        // deletions do not increase seqPos
        break;
      case(BAM_CREF_SKIP):
        // assume this is a spliced alignment for RNA, so okay
        break;
      case(BAM_CHARD_CLIP):
        // clipped region is not in seq
        break;
      case(BAM_CSOFT_CLIP):
        likelihood *= likeMiss(readQual, seqPos, count, qOff);
        seqPos += count;
        break;
    }
  }
  return likelihood;
}

// takes in a read and returns the match likelihood (due to matches, mismatches, indels)
double getMatchLikelihoodBAM(bam1_t *read, int qOff){
  assert(read != NULL);
  double likelihood;
  // read CIGAR first
  char *readQual = (char*) bam1_qual(read);
  uint32_t *cigar = bam1_cigar(read);
  int inserts = 0;
  int deletions = 0;
  int totalMatch = 0;
  //printf("getMatchLikelihoodBAM(%s, %d)\n", bam1_qname(read), qOff);

  likelihood = getCIGARLikelihoodBAM(read->core.n_cigar, cigar, readQual, qOff, &inserts, &deletions, &totalMatch);

  char *md = (char*) bam_aux_get(read, "MD");
  //printf("%s %f MD:%s\n", bam1_qname(read), likelihood, md);
  if (md != NULL && md[0] == 'Z') {
    likelihood *= getMDLikelihood(md + 1, readQual, qOff);
  } else {
    printf("WARNING: could not find the MD tag for %s\n", bam1_qname(read));
  }

  // TODO account for length of match and total length of read (soft/hard clipping).  I.e. a 10 base match is much less informative than a 20 base match
  // presently longer matches have lower likelihood because no quality is at 100%
  //printf("getMatchLikelihoodBAM(%s, %d) = %f\n", bam1_qname(read), qOff, likelihood);
  return likelihood;
}

// returns the 2-bit hash representation of a nucl. given its place in the kmer
int kmerHash(char c1, int place){
  int add1 = pow(2, 2*place);
  int add2 = add1*2;
  //printf("Checking char: %c\n", c1);
  if(c1 == 'A'){
    return 0;
  }else if(c1 == 'T'){
    return add1;
  }else if(c1 == 'C'){
    return add2;
  }else if(c1 == 'G'){
    return add1 + add2;
  }else{
    // assert(false);
    return -100000;
  }
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
  return current;
}

// apply statistics
void applyDepthAndMatchToContig(alignSet_t *alignment, contig_t *contig, double likeNormalizer) {
  double likelihood = alignment->likelihood;
  double normalLikelihood = likelihood / likeNormalizer;
  double normalLikelihood2 = likelihood; // * normalLikelihood;  // TEST!
  int j;
  if (alignment->start1 >= 0) {
    for(j = alignment->start1; j < alignment->end1; j++){
      contig->depth[j] += 1.0; //normalLikelihood; We picked a winner, it gets full prob
      contig->matchLikelihood[j] += normalLikelihood2;
    }
  }
  if (alignment->start2 >= 0) {
    for(j = alignment->start2; j < alignment->end2; j++){
      contig->depth[j] += 1.0; // normalLikelihood;
      contig->matchLikelihood[j] += normalLikelihood2;
    }
  }
}

// this applies the placement(s) to the assembly part (SINGLE PART)
// I feel like this could be sped up with a hash table vs the current linked lists, but we will see...
int applyPlacement(alignSet_t *head, assemblyT *theAssembly){

  // normalize the probs
  double likeNormalizer = getTotalLikelihood(head);

  int winner;
  alignSet_t *current = getPlacementWinner(head, likeNormalizer, &winner);

  if(current == NULL){
    //printf("No winner, failed to place %s. currentLikelihood: %f, Normalizer: %f\n", head->name, head->likelihood, likeNormalizer);
    return -1;
  }

  // apply the placement
  int i = 0;
  // this seems way too slow, why not use:
  // contig_t *contig = theAssembly->contigs[head->contigId];
  for(i = 0; i < theAssembly->numContigs; i++){ // find the right contig
    contig_t *contig = theAssembly->contigs[i];
    if(i == head->contigId){ // then add the head placement
      applyDepthAndMatchToContig(current, contig, likeNormalizer);
    }
    break;
  }
  return winner;
}

// compute the depth statistics
int computeDepthStats(assemblyT *theAssembly){
  int i, j;
  double depthNormalizer[101];
  long depthNormalizerCount[101];
  for(i = 0; i < 101; i++){
    depthNormalizer[i] = 0.0;
    depthNormalizerCount[i] = 0;
  }
  double tempLike;
  for(i = 0; i < theAssembly->numContigs; i++){ // for each contig
    contig_t *contig = theAssembly->contigs[i];
    for(j = 0; j < contig->seqLen; j++){
      float depth = contig->depth[j];
      if (depth < 0.1)
        continue;
      depthNormalizer[contig->GCcont[j]] += depth;
      depthNormalizerCount[contig->GCcont[j]] += 1;
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
      tempLike = poissonPMF(contig->depth[j], depthNormalizer[contig->GCcont[j]]);
      if(tempLike < minLogLike || isnan(tempLike)){tempLike = minLogLike;}
      contig->depthLikelihood[j] = tempLike;
      tempLike = contig->matchLikelihood[j]/contig->depth[j];
      if(tempLike < minLogLike || isnan(tempLike)){tempLike = minLogLike;}
      contig->matchLikelihood[j] = tempLike;
    }
  }
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
    int numReads = 2;
    switch(orientation) {
      case(VALID_FR):
      case(VALID_RF):
      case(VALID_FF):
        if (libParams->isSortedByName == -1) {
          libParams->isSortedByName = 1;
          printf("Setting library to be sorted by name\n");
        }
        mapLen = getMapLenBAM(thisRead);
        break;
      case(READ1_ONLY):
      case(READ2_ONLY):
        numReads = 1;
        break;
      case(UNRELATED_PAIR):
        numReads = 2;
        if (libParams->isSortedByName == -1) {
          libParams->isSortedByName = 0;
          printf("Setting library to be unsorted by name\n");
        }
        break;
      default:
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
    printf("Evaluating %s orientation\n", MATE_ORIENTATION_LABELS[j]);
    libraryMateParametersT *mateParams = &libParams->mateParameters[j];
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
    if (j == READ1_ONLY || j == READ2_ONLY)
      mateParams->libraryFraction /= 2.0;

    if (j == VALID_FR || j == VALID_RF || j == VALID_FF || j == SINGLE_READ) {

      // TODO better test for significance and normal distribution for a valid orientation
      if (mateParams->libraryFraction > 0.02) {
        totalValid += mateParams->count * (j == SINGLE_READ ? 1.0 : 2.0);
        mateParams->isValid = 1;
      }

      printf("Read %ld properly oriented sequences (%s) and purged %ld %0.1lf%% & %0.1lf%% outliers\n",
          mateParams->count, MATE_ORIENTATION_LABELS[j], purged, (outlierFraction*100.0), ((1.0-outlierFraction)*100.0));
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
  printf("There were %ld improper pairs\n", improperPair);
  libParams->avgReadSize = libParams->avgReadSize / libParams->numReads;
  printf("Found sample avg read size to be %ld\n", libParams->avgReadSize);

  double totalSingle = libParams->mateParameters[READ1_ONLY].count + libParams->mateParameters[READ2_ONLY].count + libParams->mateParameters[SINGLE_READ].count + libParams->mateParameters[UNMAPPED_SINGLE].count;
  libParams->totalSingleFraction = totalSingle / libParams->numReads;
  libParams->totalValidFraction = totalValid / libParams->numReads;
  libParams->totalChimerFraction = 1.0 - libParams->totalSingleFraction - libParams->totalValidFraction;

  printf("Valid %0.3lf%%, Single %0.3lf%%, Chimer %0.3lf%%\n", libParams->totalValidFraction*100, libParams->totalSingleFraction*100, libParams->totalChimerFraction*100);

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
  read2Only->start1 = thisAlignment->start2;
  read2Only->end1 = thisAlignment->end2;

  // reset both alignments read2 coordinates
  thisAlignment->start2 = thisAlignment->end2 = -1;
  read2Only->start2 = read2Only->end2 = -1;

  read2Only->likelihood = likelihood;
  strncpy(read2Only->name, bam1_qname(thisReadMate), MAX_NAME_LENGTH);
  read2Only->name[MAX_NAME_LENGTH] = '\0'; // ensure null termination
  if (read2Only->likelihood > 0.0) {
    read2Only->contigId = thisReadMate->core.tid;
  }
  read2Only->nextAlignment = NULL;
}

enum MATE_ORIENTATION setAlignment(bam_header_t *header, assemblyT *theAssembly, alignSet_t *thisAlignment, alignSet_t *secondaryAlignment, void **mateTree, libraryParametersT *libParams, enum MATE_ORIENTATION orientation, bam1_t *thisRead, bam1_t *thisReadMate) {
  assert(thisAlignment != NULL);
  double likelihoodRead1 = 1.0;
  double likelihoodRead2 = 1.0;
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
    likelihoodRead1  = getMatchLikelihoodBAM(thisRead, qOff);
    thisAlignment->start1 = thisRead->core.pos;
    thisAlignment->end1   = bam_calend(&thisRead->core, bam1_cigar(thisRead));
    assert(thisAlignment->start1 <= thisAlignment->end1);
    thisAlignment->contigId = thisRead->core.tid;
  }

  if (thisReadMate == NULL || (thisReadMate->core.flag & BAM_FUNMAP) == BAM_FUNMAP) {
    thisAlignment->start2 = -1;
    thisAlignment->end2   = -1;
  } else {
    likelihoodRead2  = getMatchLikelihoodBAM(thisReadMate, qOff);
    thisAlignment->start2 = thisReadMate->core.pos;
    thisAlignment->end2   =  bam_calend(&thisReadMate->core, bam1_cigar(thisReadMate));
    assert(thisAlignment->start2 <= thisAlignment->end2);
    if (thisAlignment->start1 < 0) {
      thisAlignment->contigId = thisReadMate->core.tid;
    }
  }

  libraryMateParametersT *mateParameters = &libParams->mateParameters[orientation];

  switch (orientation) {
    case (VALID_FR):
    case (VALID_RF):
    case (VALID_FF):
      // two reads

      if (mateParameters->isValid) {
        // valid orientation
        likelihoodInsert = getInsertLikelihoodBAM(thisRead, mateParameters->insertLength, mateParameters->insertStd);
        thisAlignment->likelihood = libParams->totalValidFraction * likelihoodInsert * likelihoodRead1 * likelihoodRead2;
        break;
      } else {
        // change the orientation... this is actually a chimer
        if (thisRead->core.tid == thisRead->core.mtid)
          orientation = CHIMER_SAME_CONTIG;
        else
          orientation = CHIMER_DIFF_CONTIG;
      }
      //  continue... this is actually a chimer
    case (CHIMER_SAME_CONTIG) :
    case (CHIMER_DIFF_CONTIG) :
      //printf("WARNING: chimeric read mate pair %s.\n", bam1_qname(thisRead));


      likelihoodInsert = libParams->totalChimerFraction;
      if (orientation == CHIMER_DIFF_CONTIG) {
        // TODO refine based on proximity to end of contigs
        // i.e. factor in: likelihoodThatRead1AndRead2AreCloseToContigEdgeSoAreNotChimersButProbablySpanningMatePairs
      }

      thisAlignment->likelihood *= likelihoodInsert * likelihoodRead1;
      double likelihoodMate = likelihoodInsert * likelihoodRead2;

      // set secondaryAlignment to map both reads separately...
      if (thisReadMate != NULL)
        setSingleRead2Alignment(header, secondaryAlignment, thisAlignment, thisReadMate, likelihoodMate);

      break;

    case (UNRELATED_PAIR):
    case (READ1_ONLY):
    case (READ2_ONLY):

      if(thisReadMate != NULL && (thisReadMate->core.flag & BAM_FUNMAP) != BAM_FUNMAP) {
        // set seconaryAlignment to thisReadMate, remove thisReadMate from thisAlignment

        setSingleRead2Alignment(header, secondaryAlignment, thisAlignment, thisReadMate, likelihoodRead2);

        if ((thisReadMate->core.flag & BAM_FMUNMAP) == BAM_FMUNMAP) {
          // mate is not mapped, apply likelihood now
          secondaryAlignment->likelihood *= libParams->totalSingleFraction;
        } else {
          assert(thisReadMate->core.tid == thisReadMate->core.mtid);
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
        thisAlignment->likelihood = likelihoodRead1;
        thisAlignment->start2 = thisAlignment->end2 = -1;
        if ((thisRead->core.flag & BAM_FMUNMAP) == BAM_FMUNMAP) {
          // mate is not mapped, apply likelihood now
          likelihoodInsert = libParams->totalSingleFraction;
          thisAlignment->likelihood *= likelihoodInsert;
        } else {
          assert(thisRead->core.tid == thisRead->core.mtid);
          // store this or get mate if already seen
          alignSet_t *mateAlignment = getOrStoreMateAlignment(mateTree, thisAlignment, thisRead);
          if (mateAlignment != NULL) {
            likelihoodRead2 = mateAlignment->likelihood;
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
      thisAlignment->likelihood = libParams->totalSingleFraction * likelihoodRead1;
      if (thisReadMate != NULL)
        setSingleRead2Alignment(header, secondaryAlignment, thisAlignment, thisReadMate, libParams->totalSingleFraction * likelihoodRead2);
      break;

    case (UNMAPPED_PAIR):
      break;

    default :
      assert(thisRead != NULL && thisReadMate != NULL);
      printf("Skipping %s read %s %s\n", MATE_ORIENTATION_LABELS[orientation], bam1_qname(thisRead), bam1_qname(thisReadMate));
      thisAlignment->likelihood = 0.0;
      break;
  }

  return orientation;
}

// NORMALIZATION
// divide by the expected likelihood of the read by the normalization factor Z (from Bayes rule)
// given only its length and the parameters of the distributions (See paper appendix)
double zNormalization(libraryMateParametersT *mateParams, bam1_t *thisRead, bam1_t *thisReadMate, int qOff){
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

  // find the expected match score
  double expMatch = 0.0;
  for(i = 0; i < totalLen; i++){
    expMatch += pow(pow(Qavg, i)*pow((1.0 - Qavg)/3.0, totalLen - i - 1), 2);
  }

  //printf("Qavg: %f\n", Qavg);
  //printf("expMatch: %f\n", expMatch);

  // int((pmf of normal(0,sigma))^2, 0, infty)
  double expIns = 1.0/(2*sqrt(3.14159265)*mateParams->insertStd);

  return expMatch*expIns;
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
    // given only its length and the parameters of the distributions (See paper appendix)
    thisAlignment->likelihood /= zNormalization(&libParams->mateParameters[orientation], thisRead, libParams->isSortedByName == 1 ? thisReadMate : NULL, libParams->qOff);
    

    if (orientation == NO_READS)
      break;

    //printf("%s\n", bam1_qname(thisRead));

    //process secondaryAlignment (thisReadMate separately)
    // apply placement of read2 to target2...
    // HACK!!!
    // TODO fix for multiple possible placements
    if (thisReadMate != NULL && secondaryAlignment.likelihood > 0.0) {
      int winner = applyPlacement(&secondaryAlignment, theAssembly);
      if (winner < 0) {
        printf("WARNING: no placement found for read2 of chimer %s!\n", secondaryAlignment.name);
      } else {
        if (placementBam != NULL)
          bam_write1(placementBam->x.bam, thisReadMate);
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
