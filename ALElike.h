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

const static double lnfactconst2 = 0.918938533204672741780329;
const static double minLogLike = -60.0;

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

double getInsertLikelihood(SAM_t *read, double mu, double var){
    double likelihood = GetInsertProbNormal(abs((double)read->mapLen) - mu, var);
    //printf("Insert(%f): %f drawn from N(%f,%f)? = %lf\n", (double)read->mapLen, abs((double)read->mapLen) - mu, 0.0, var, likelihood);
    return likelihood;
}

double getInsertLikelihoodBAM(bam1_t *read1, bam1_t *read2, double mu, double var){
	double likelihood = GetInsertProbNormal(abs(getMapLenBAM(read1,read2) - mu), var);
	return likelihood;
}

// finds the likelihood of a string of misses in read from seqPos to matchLen
double likeMiss(char *readQual, int seqPos, int missLen, int qOff){
    int i;
    double likelihood = 1.0;
    for(i = seqPos; i < seqPos + missLen; i++){
        //likelihood = likelihood*((1.0 - 0.99)/3.0);
        likelihood = likelihood*((1.0 - QtoP[readQual[i] - qOff])/3.0); // sometimes want -33/64?
    }
    return likelihood;
}

// finds the likelihood of a string of matches in read from seqPos to matchLen
double likeMatch(char *readQual, int seqPos, int matchLen, int qOff){
    int i;
    double likelihood = 1.0;
    for(i = seqPos; i < seqPos + matchLen; i++){
        //likelihood = likelihood*(0.99);
        likelihood = likelihood*(QtoP[readQual[i] - qOff]);// sometimes want -33/64?
    	//printf("likeMatch %d %c %f %f\n", i, readQual[i], QtoP[readQual[i] - qOff], likelihood);
    }
    //printf("likeMatch: %s %lf\n", read->readName, likelihood);
    return likelihood;
}

double getMDLikelihood(char *MD, char *readQual, int qOff) {
	int stop = 0;
	int pos = 0;
	int seqPos = 0;
	int matchLen, missLen, delLen;
	int totalMatch = 0, totalMiss = 0, totalDel = 0;
	double likelihood = 1.0;
	//printf("getMDLikelihood(%s, %s, %d)\n", MD, readQual, qOff);
	// parse MD field
	while(stop == 0){
		// matches
		matchLen = 0;
		while(isdigit(MD[pos])){
			matchLen = matchLen*10 + hackedIntCast(MD[pos]);
			pos++;
		}
		if (matchLen > 0) {
			totalMatch += matchLen;
			likelihood = likelihood*likeMatch(readQual, seqPos, matchLen, qOff);
			//printf("MD %d match %d. %f\n", seqPos, matchLen, likelihood);
			seqPos += matchLen;
		}
		// misses
		missLen = 0;
		while(MD[pos] == 'A' || MD[pos] == 'T' || MD[pos] == 'C' || MD[pos] == 'G' || MD[pos] == '0'){
			if(MD[pos] == 'A' || MD[pos] == 'T' || MD[pos] == 'C' || MD[pos] == 'G'){
				missLen++;
				likelihood = likelihood*likeMiss(readQual, seqPos, 1, qOff);
				seqPos++;
			}
			if(MD[pos] == 'N'){
				missLen++;
				likelihood = likelihood*0.25;
				seqPos++;
			}
			pos++;
			//printf("MD %d miss  %d. %f\n", seqPos, 1, likelihood);
		}
		totalMiss += missLen;
		// deletions
		delLen = 0;
		if(MD[pos] == '^'){
			pos++;
			while(isalpha(MD[pos])){
				delLen++;
				pos++;
				// assume likelihood of deletion is same as substitution
				// TODO revise to be targeting proper base quality score (or avg over two)
				likelihood = likelihood*likeMiss(readQual, seqPos, 1, qOff);
				//printf("MD %d del   %d. %f\n", seqPos, delLen, likelihood);
			}
		}
		totalDel += delLen;

		// sees if we are at the end
		if(MD[pos] == '\0'){
			stop = 1;
		}
	}
	if (totalMatch == 0)
		likelihood = 0.0;

	return likelihood;

}
// takes in a read and returns the match likelihood (due to matches, mismatches, indels)
double getMatchLikelihood(SAM_t *read, int qOff){

    int stop = 0;
    int pos = 5;
    int totalDel = 0, totalIns = 0;

    // parse MD
    double likelihood = getMDLikelihood(&read->MD[5], read->readQual, qOff);

     // parse CIGAR
    stop = 0;
    pos = 0;
    int num;
    while(stop == 0){
        // matches
        num = 0;
        while(isdigit(read->cigar[pos])){
            num = num*10 + hackedIntCast(read->cigar[pos]);
            pos++;
        }
        // insertions
        if(read->cigar[pos] == 'I'){
        	// TODO assume likelihood of insertion is same as substitution
            totalIns += num;
        }
        pos++;
        // read to the end
        if(read->cigar[pos] == '\0'){
            stop = 1;
        }
    }
    if(totalDel + totalIns > 0){
        likelihood = 0.0; // can modify later to include insertions and deletions
    }
    return likelihood;
}

// takes in a read and returns the match likelihood (due to matches, mismatches, indels)
double getMatchLikelihoodBAM(bam1_t *read, int qOff){
	double likelihood = 1.0;
	// read CIGAR first
	int inserts = 0;
	int deletions = 0;
	int totalMatch = 0;
	uint32_t *cigar = bam1_cigar(read);
	int i;
	for(i=0 ; i < read->core.n_cigar ; i++) {
		int32_t count = (*cigar >> BAM_CIGAR_SHIFT);
		switch (*cigar & BAM_CIGAR_MASK) {
		case(BAM_CMATCH) : totalMatch += count; break;;
		case(BAM_CINS)   : inserts += count; break;;
		case(BAM_CDEL)    : deletions += count; break;;
		}
	}
	if (deletions+inserts > 0) {
		printf("Detected indel on %s\n", bam1_qname(read));
		likelihood = 0.0;  // can modify later to include insertions and deletions
	} else {
		char *md = (char*) bam_aux_get(read, "MD");
		//printf("%s %f MD:%s\n", bam1_qname(read), likelihood, md);
		if (md != NULL && md[0] == 'Z') {
			likelihood *= getMDLikelihood(md + 1, (char*) bam1_qual(read), qOff);
		} else {
			printf("WARNING: could not find the MD tag for %s\n", bam1_qname(read));
		}
	}
	return likelihood;
}

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
        return -100000;
    }
}

int getKmerHash(char *seq, int startPos, int kmerLen){
    int i;
    int hash = 0;
    for(i = 0; i < kmerLen; i++){
        hash += kmerHash(seq[startPos+i], i);
    }
    return hash;
}

void computeKmerStats(assemblyT *theAssembly, int kmer){
    int i, j, k, totalKmers, hash;
    // calculate total possible kmers
    totalKmers = 1;
    for(i = 0; i < kmer; i++){
        totalKmers = totalKmers*4;
    }
    int *kmerVec = malloc(totalKmers*sizeof(int));
    // find all kmers present
    if(theAssembly->numContigs > 1){
        for(i = 0; i < theAssembly->numContigs; i++){
            // initialize kmerVec
            for(j = 0; j < totalKmers; j++){
                kmerVec[j] = 0;
            }
            totalKmers = 0;
            // add up the kmers
            for(j = 0; j < theAssembly->contigs[i].seqLen - kmer + 1; j++){
                hash = getKmerHash(theAssembly->contigs[i].seq, j, kmer);
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
                theAssembly->contigs[i].kmerLikelihood[j] = 0.0;
                for(k = 0; k < j+1; k++){
                    hash = getKmerHash(theAssembly->contigs[i].seq, k, kmer);
                    if(hash > -1){
                        theAssembly->contigs[i].kmerLikelihood[j] = theAssembly->contigs[i].kmerLikelihood[j] + 1.0/(float)(j+1)*(float)(kmerVec[hash])/(float)(totalKmers);
                    }
                }
                //printf("New likelihood[%i]: %f.\n", j, theAssembly->contigs[i].kmerLikelihood[j]);
            }
            //printf("First.\n");
            // middle bunch
            for(j = kmer; j < theAssembly->contigs[i].seqLen - kmer; j++){
                for(k = 0; k < kmer; k++){
                    hash = getKmerHash(theAssembly->contigs[i].seq, j - k, kmer);
                    if(hash > -1){
                        theAssembly->contigs[i].kmerLikelihood[j] = theAssembly->contigs[i].kmerLikelihood[j] + 1.0/(float)(kmer)*(float)(kmerVec[hash])/(float)(totalKmers);
                    }
                }
                //printf("New likelihood[%i]: %f.\n", j, theAssembly->contigs[i].kmerLikelihood[j]);
            }
            //printf("Mid.\n");
            // last bits
            for(j = theAssembly->contigs[i].seqLen - kmer; j < theAssembly->contigs[i].seqLen; j++){
                theAssembly->contigs[i].kmerLikelihood[j] = 0.0;
                for(k = j - kmer + 1; k < j - kmer + 1 + (theAssembly->contigs[i].seqLen - j); k++){
                    hash = getKmerHash(theAssembly->contigs[i].seq, k, kmer);
                    if(hash > -1){
                        theAssembly->contigs[i].kmerLikelihood[j] = theAssembly->contigs[i].kmerLikelihood[j] + 1.0/(float)(theAssembly->contigs[i].seqLen - j)*(float)(kmerVec[hash])/(float)(totalKmers);
                    }
                }
                //printf("New likelihood[%i]: %f.\n", j, theAssembly->contigs[i].kmerLikelihood[j]);
            }
            //printf("Last.\n");
            totalKmers = 0;
        }
    }else{ // ONLY ONE CONTIG IN THE ASSEMBLY!!
        // initialize kmerVec
        for(j = 0; j < totalKmers; j++){
            kmerVec[j] = 0;
        }
        totalKmers = 0;
        // add up the kmers
        for(j = 0; j < theAssembly->contigs->seqLen - kmer + 1; j++){
            hash = getKmerHash(theAssembly->contigs->seq, j, kmer);
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
            theAssembly->contigs->kmerLikelihood[j] = 0.0;
            for(k = 0; k < j+1; k++){
                hash = getKmerHash(theAssembly->contigs->seq, k, kmer);
                if(hash > -1){
                    theAssembly->contigs->kmerLikelihood[j] = theAssembly->contigs->kmerLikelihood[j] + 1.0/(float)(j+1)*(float)(kmerVec[hash])/(float)(totalKmers);
                }
            }
            //printf("New likelihood[%i]: %f.\n", j, theAssembly->contigs[i].kmerLikelihood[j]);
        }
        //printf("First.\n");
        // middle bunch
        for(j = kmer; j < theAssembly->contigs->seqLen - kmer; j++){
            for(k = 0; k < kmer; k++){
                hash = getKmerHash(theAssembly->contigs->seq, j - k, kmer);
                //printf("Hash: %i\n", hash);
                if(hash > -1){
                    theAssembly->contigs->kmerLikelihood[j] = theAssembly->contigs->kmerLikelihood[j] + 1.0/(float)(kmer)*(float)(kmerVec[hash])/(float)(totalKmers);
                    //printf("1.0/%f*%f/%f = %f\n", (float)kmer, (float)kmerVec[hash], (float)totalKmers, 1.0/(float)(kmer)*(float)(kmerVec[hash])/(float)(totalKmers));
                }
            }
            //printf("New likelihood[%i]: %f.\n", j, theAssembly->contigs->kmerLikelihood[j]);
        }
        //printf("Mid.\n");
        // last bits
        for(j = theAssembly->contigs->seqLen - kmer; j < theAssembly->contigs->seqLen; j++){
            theAssembly->contigs->kmerLikelihood[j] = 0.0;
            for(k = j - kmer + 1; k < j - kmer + 1 + (theAssembly->contigs->seqLen - j); k++){
                hash = getKmerHash(theAssembly->contigs->seq, k, kmer);
                if(hash > -1){
                    theAssembly->contigs->kmerLikelihood[j] = theAssembly->contigs->kmerLikelihood[j] + 1.0/(float)(theAssembly->contigs->seqLen - j)*(float)(kmerVec[hash])/(float)(totalKmers);
                }
            }
            //printf("New likelihood[%i]: %f.\n", j, theAssembly->contigs[i].kmerLikelihood[j]);
        }
        //printf("Last.\n");
        totalKmers = 0;
    }
}

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

alignSet_t *getPlacementWinner(alignSet_t *head, double likeNormalizer, int *winner) {

	*winner = -1;
	if(likeNormalizer == 0.0){ // no real placement
		return NULL;
	}

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

void applyDepthAndMatchToAssembly(alignSet_t *alignment, assemblyT *theAssembly, int contigId, double likeNormalizer) {
	contig_t *contig = &(theAssembly->contigs[contigId]);
	double likelihood = alignment->likelihood;
	double normalLikelihood = likelihood / likeNormalizer;
	double normalLikelihood2 = likelihood * normalLikelihood;
	int j;
	if (alignment->start1 >= 0) {
		for(j = alignment->start1; j < alignment->end1; j++){
			contig->depth[j] += normalLikelihood;
			contig->matchLikelihood[j] += normalLikelihood2;
		}
	}
	if (alignment->start2 >= 0) {
		for(j = alignment->start2; j < alignment->end2; j++){
			contig->depth[j] += normalLikelihood;
			contig->matchLikelihood[j] += normalLikelihood2;
		}
	}
}
// this applies the placement(s) to the assembly part (SINGLE PART)
// I feel like this could be sped up with a hash table vs the current linked lists, but we will see...
int applyPlacement(alignSet_t *head, assemblyT *theAssembly){

	// instead of random choice, choose a consistant, but random choice
	unsigned int iseed = JSHash(head->name);
	srand(iseed); //srand ((unsigned int) time(NULL));

	// normalize the probs
	double likeNormalizer = getTotalLikelihood(head);

	int winner;
	alignSet_t *current = getPlacementWinner(head, likeNormalizer, &winner);

	if(current == NULL){
		printf("No winner, failed to place %s. currentLikelihood: %f, Normalizer: %f\n", head->name, head->likelihood, likeNormalizer);
		return -1;
	}

	// apply the placement
	int i = 0;
	for(i = 0; i < theAssembly->numContigs; i++){ // find the right contig
		if(strcmp(theAssembly->contigs[i].name, head->mapName) == 0){ // then add the head placement
			applyDepthAndMatchToAssembly(current, theAssembly, i, likeNormalizer);
		}
		break;
	}
	return winner;
}

// this applies the placement(s) to the assembly part(s)
// I feel like this could be sped up with a hash table vs the current linked lists, but we will see...
int applyPlacementOld(alignSet_t *head, assemblyT *theAssembly){
    // normalize the probs
    double likeNormalizer = 0.0;
    likeNormalizer += head->likelihood;
    alignSet_t *current = head;
    while(current->nextAlignment != NULL){
        current = current->nextAlignment;
        likeNormalizer += current->likelihood;
    }
    //printf("Normalizer: %f\n", likeNormalizer);
    if(likeNormalizer == 0.0){ // no real placement
	//printf("Failed to place\n");
        return -1;
    }
    
    // apply the first placement
    //contig_t *matchContig;
    int i, j;
    if(theAssembly->numContigs > 1){
        for(i = 0; i < theAssembly->numContigs; i++){ // find the right contig
            if(strcmp(theAssembly->contigs[i].name, head->mapName) == 0){ // then add the head placement
                for(j = head->start1; j < head->end1; j++){
                    theAssembly->contigs[i].depth[j] = theAssembly->contigs[i].depth[j] + head->likelihood/likeNormalizer;
                    theAssembly->contigs[i].matchLikelihood[j] += head->likelihood*(head->likelihood/likeNormalizer);
                }
                for(j = head->start2; j < head->end2; j++){
                    theAssembly->contigs[i].depth[j] = theAssembly->contigs[i].depth[j] + head->likelihood/likeNormalizer;
                    theAssembly->contigs[i].matchLikelihood[j] += head->likelihood*(head->likelihood/likeNormalizer);
                }
                //break;
            }
        }
        // do the rest
        current = head;
        while(current->nextAlignment != NULL){
            current = current->nextAlignment;
            for(i = 0; i < theAssembly->numContigs; i++){ // find the right contig
                if(strcmp(theAssembly->contigs[i].name, current->mapName) == 0){ // then add the head placement
                    for(j = current->start1; j < current->end1; j++){
                        theAssembly->contigs[i].depth[j] = theAssembly->contigs[i].depth[j] + current->likelihood/likeNormalizer;
                        theAssembly->contigs[i].matchLikelihood[j] += current->likelihood*(current->likelihood/likeNormalizer);
                    }
                    for(j = current->start2; j < current->end2; j++){
                        theAssembly->contigs[i].depth[j] = theAssembly->contigs[i].depth[j] + current->likelihood/likeNormalizer;
                        theAssembly->contigs[i].matchLikelihood[j] += current->likelihood*(current->likelihood/likeNormalizer);
                    }
                    //break;
                }
            }
        }
    }else{
        if(strcmp(theAssembly->contigs->name, head->mapName) == 0){ // then add the head placement
            for(j = head->start1; j < head->end1; j++){
                theAssembly->contigs->depth[j] = theAssembly->contigs->depth[j] + head->likelihood/likeNormalizer;
                theAssembly->contigs->matchLikelihood[j] += head->likelihood*(head->likelihood/likeNormalizer);
		//if(j == 3528300){printf("depthf %i = %lf, like = %lf, norm = %lf\n", j, theAssembly->contigs->depth[j], head->likelihood, likeNormalizer);printAlignments(head);}
            }
            for(j = head->start2; j < head->end2; j++){
                theAssembly->contigs->depth[j] = theAssembly->contigs->depth[j] + head->likelihood/likeNormalizer;
                theAssembly->contigs->matchLikelihood[j] += head->likelihood*(head->likelihood/likeNormalizer);
		//if(j == 3528300){printf("deptht %i = %lf, like = %lf, norm = %lf\n", j, theAssembly->contigs->depth[j], head->likelihood, likeNormalizer);printAlignments(head);}
            }
        }
        // do the rest
        current = head;
        while(current->nextAlignment != NULL){
            current = current->nextAlignment;
            if(strcmp(theAssembly->contigs->name, current->mapName) == 0){ // then add the head placement
                for(j = current->start1; j < current->end1; j++){
                    theAssembly->contigs->depth[j] = theAssembly->contigs->depth[j] + current->likelihood/likeNormalizer;
                    theAssembly->contigs->matchLikelihood[j] += current->likelihood*(current->likelihood/likeNormalizer);
		    //if(j == 3528300){printf("depth2f %i = %lf, like = %lf, norm = %lf\n",j,  theAssembly->contigs->depth[j], current->likelihood, likeNormalizer);printAlignments(head);}
                }
                for(j = current->start2; j < current->end2; j++){
                    theAssembly->contigs->depth[j] = theAssembly->contigs->depth[j] + current->likelihood/likeNormalizer;
                    theAssembly->contigs->matchLikelihood[j] += current->likelihood*(current->likelihood/likeNormalizer);
		    //if(j == 3528300){printf("depth2t %i = %lf, like = %lf, norm = %lf\n", j, theAssembly->contigs->depth[j], current->likelihood, likeNormalizer);printAlignments(head);}
                }
            }
        }
    }
    return 0;
}

int computeDepthStats(assemblyT *theAssembly){
    int i, j;
    double depthNormalizer[100];
    int depthNormalizerCount[100];
    for(i = 0; i < 100; i++){
	depthNormalizer[i] = 0.0;
	depthNormalizerCount[i] = 0;
    }
    double tempLike;
    if(theAssembly->numContigs > 1){
        for(i = 0; i < theAssembly->numContigs; i++){ // for each contig
            for(j = 0; j < theAssembly->contigs[i].seqLen; j++){
		//printf("%f %i\n", 100.0*theAssembly->contigs[i].GCcont[j], (int)floor(100.0*theAssembly->contigs[i].GCcont[j]));
                depthNormalizer[(int)floor(100.0*theAssembly->contigs[i].GCcont[j])] += theAssembly->contigs[i].depth[j];
		depthNormalizerCount[(int)floor(100.0*theAssembly->contigs[i].GCcont[j])] += 1;
            }
            for(j = 0; j < 100; j++){
	      if(depthNormalizerCount[j] > 0){
		depthNormalizer[j] = depthNormalizer[j]/(double)depthNormalizerCount[j];
	      }else{
		depthNormalizer[j] = 0.1;
	      }
	    }
            for(j = 0; j < theAssembly->contigs[i].seqLen; j++){
                tempLike = poissonPMF(theAssembly->contigs[i].depth[j], depthNormalizer[(int)floor(100.0*theAssembly->contigs[i].GCcont[j])]);
                if(tempLike < minLogLike || isnan(tempLike)){tempLike = minLogLike;}
                theAssembly->contigs[i].depthLikelihood[j] = tempLike;
                theAssembly->contigs[i].matchLikelihood[j] = theAssembly->contigs[i].matchLikelihood[j]/theAssembly->contigs[i].depth[j];
            }
        }
    }else{
        for(j = 0; j < theAssembly->contigs->seqLen; j++){
	    //printf("%f %i\n", 100.0*theAssembly->contigs->GCcont[j], (int)floor(100.0*theAssembly->contigs->GCcont[j]));
                depthNormalizer[(int)floor(100.0*theAssembly->contigs->GCcont[j])] += theAssembly->contigs->depth[j];
		depthNormalizerCount[(int)floor(100.0*theAssembly->contigs->GCcont[j])] += 1;
            }
        for(j = 0; j < 100; j++){
	      if(depthNormalizerCount[j] > 0){
		depthNormalizer[j] = depthNormalizer[j]/(double)depthNormalizerCount[j];
	      }else{
		depthNormalizer[j] = 0.1;
	      }
	    }
        for(j = 0; j < theAssembly->contigs->seqLen; j++){
            tempLike = poissonPMF(theAssembly->contigs->depth[j], depthNormalizer[(int)floor(100.0*theAssembly->contigs->GCcont[j])]);
            if(tempLike < minLogLike || isnan(tempLike)){tempLike = minLogLike;}
            theAssembly->contigs->depthLikelihood[j] = tempLike;
            theAssembly->contigs->matchLikelihood[j] = theAssembly->contigs->matchLikelihood[j]/theAssembly->contigs->depth[j];
        }
    }
    return 1;
}
