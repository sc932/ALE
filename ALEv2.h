struct contig_struct{
    char name[256];
    int seqLen;
    char *seq;
    double *depth;
    double *matchLikelihood;
    double *kmerLikelihood;
    double *depthLikelihood;
};

struct assembly_struct{
    int numContigs;
    int totalAssemLen;
    struct contig_struct *contigs;
};

struct SAMalignment{ // values in order of read in
    char readName[256];
    int outInfo;
    char refName[256];
    int mapStart;
    int mapPair;
    char cigar[256];
    char flag2[10];
    int mapEnd;
    int mapLen;
    char readSeq[256];
    char readQual[256];
    char XA[256];
    char MD[256];
    char NM[256];
    struct SAMalignment *pair;
};

struct setOfAlignments{
    double likelihood;
    int start1, start2;
    int end1, end2;
    char name[256];
    char mapName[256];
    struct setOfAlignments *nextAlignment;
};

typedef struct setOfAlignments alignSet_t;
typedef struct SAMalignment SAM_t;
typedef struct contig_struct contig_t;
typedef struct assembly_struct assemblyT;

void printAssembly(assemblyT *theAssembly);

// int sizeOfContig(int len){
//     return 256*sizeof(char) + sizeof(int) + len*(sizeof(char) + 5*sizeof(double));
// }

int findNumAssemPieces(kseq_t *ins){
    int l, count = 0;
    while ((l = kseq_read(ins)) >= 0) {
        count++;
    }
    return count;
}

void readAssembly(kseq_t *ins, assemblyT *theAssembly){
    int contigLen, i, l, j = 0;
    while ((l = kseq_read(ins)) >= 0) {
        contigLen = (int)(ins->seq.l);
        printf("Found contigLen = %i\n", contigLen);
        if(theAssembly->numContigs > 1){
            theAssembly->contigs[j].seqLen = contigLen;
            theAssembly->contigs[j].seq = malloc(contigLen*sizeof(char));
            theAssembly->contigs[j].depth = malloc(contigLen*sizeof(double));
            theAssembly->contigs[j].matchLikelihood = malloc(contigLen*sizeof(double));
            theAssembly->contigs[j].kmerLikelihood = malloc(contigLen*sizeof(double));
            theAssembly->contigs[j].depthLikelihood = malloc(contigLen*sizeof(double));
            strcpy(theAssembly->contigs[j].name, ins->name.s);
            for(i = 0; i < contigLen; i++){
                theAssembly->contigs[j].seq[i] = toupper(ins->seq.s[i]);
                theAssembly->contigs[j].depth[i] = 0.0;
                theAssembly->contigs[j].matchLikelihood[i] = 0.0;
                theAssembly->contigs[j].kmerLikelihood[i] = 0.0;
                theAssembly->contigs[j].depthLikelihood[i] = 0.0;
            }
        }else{
            theAssembly->contigs->seqLen = contigLen;
            theAssembly->contigs->seq = malloc(contigLen*sizeof(char));
            theAssembly->contigs->depth = malloc(contigLen*sizeof(double));
            theAssembly->contigs->matchLikelihood = malloc(contigLen*sizeof(double));
            theAssembly->contigs->kmerLikelihood = malloc(contigLen*sizeof(double));
            theAssembly->contigs->depthLikelihood = malloc(contigLen*sizeof(double));
            strcpy(theAssembly->contigs->name, ins->name.s);
            for(i = 0; i < contigLen; i++){
                theAssembly->contigs->seq[i] = toupper(ins->seq.s[i]);
                theAssembly->contigs->depth[i] = 0.0;
                theAssembly->contigs->matchLikelihood[i] = 0.0;
                theAssembly->contigs->kmerLikelihood[i] = 0.0;
                theAssembly->contigs->depthLikelihood[i] = 0.0;
            }
        }
        j++;
    }
}

void printAssembly(assemblyT *theAssembly){
    int i;
    if(theAssembly->numContigs > 1){
        for(i = 0; i < theAssembly->numContigs; i++){
            printf("Contig %i: %s: %s\n",i, theAssembly->contigs[i].name, theAssembly->contigs[i].seq);
        }
    }else{
        printf("Contig 0: %s: %s\n", theAssembly->contigs->name, theAssembly->contigs->seq);
    }
}

// prints out the SAM alignment
int printSAM(SAM_t read){
    printf("SAM alignment:\n");
    printf("%s %i %s %i %i %s %s %i %i %s %s %s %s %s\n", read.readName, read.outInfo, read.refName, read.mapStart, read.mapPair, read.cigar, read.flag2, read.mapEnd, read.mapLen, read.readSeq, read.readQual, read.XA, read.MD, read.NM);
}

// returns the length of a sequence, does not take indels into account
int getSeqLen(char seq[256]){
//     int len = 0;
//     while(seq[len] != '\0'){
//         len++;
//     }
//     return len;
    int i;
    for(i = 0; i < 256; i++){
        if(seq[i] != 'A' && seq[i] != 'T' && seq[i] != 'C' && seq[i] != 'G'){
            //printf("Len: %i\n", i);
            return i;
        }
    }
    printf("wtf?");
    return 0;
}

// prints out all of the alignments in the linked list
void printAlignments(alignSet_t *head){
    // print out the head
    printf("Alignment 1 for read %s: %f at %i-%i and %i-%i.\n", head->name, head->likelihood, head->start1, head->end1, head->start2, head->end2);
    alignSet_t *current = head;
    int i = 1;
    while(current->nextAlignment != NULL){
        current = current->nextAlignment;
        i++;
        printf("Alignment %i for read %s: %f at %i-%i and %i-%i.\n", i, current->name, current->likelihood, current->start1, current->end1, current->start2, current->end2);
    }
}

void writeToOutput(assemblyT *theAssembly, FILE *out){
    int i, j;
    if(theAssembly->numContigs > 1){
        for(i = 0; i < theAssembly->numContigs; i++){
            fprintf(out, ">%s : depth : ln(depthLike) : placeLike : kmerLike : length=%i\n", theAssembly->contigs[i].name, theAssembly->contigs[i].seqLen);
            for(j = 0; j < theAssembly->contigs[i].seqLen; j++){
                fprintf(out, "%f,%f,%f,%f\n", theAssembly->contigs[i].depth[j], theAssembly->contigs[i].depthLikelihood[j], theAssembly->contigs[i].matchLikelihood[j], theAssembly->contigs[i].kmerLikelihood[j]);
            }
        }
    }else{
        fprintf(out, ">%s %i > depth ln(depthLike) ln(placeLike) ln(kmerLike) ln(totalLike)\n", theAssembly->contigs->name, theAssembly->contigs->seqLen);
        for(j = 0; j < theAssembly->contigs->seqLen; j++){
            fprintf(out, "%lf %lf %lf %lf %lf\n", theAssembly->contigs->depth[j], theAssembly->contigs->depthLikelihood[j], log(theAssembly->contigs->matchLikelihood[j]), log(theAssembly->contigs->kmerLikelihood[j]), theAssembly->contigs->depthLikelihood[j] + log(theAssembly->contigs->matchLikelihood[j]) + log(theAssembly->contigs->kmerLikelihood[j]));
        }
    }
}

int assemblySanityCheck(assemblyT *theAssembly){
    int i, num = theAssembly->numContigs;
    if(num == 1){
        for(i = 0; i < theAssembly->contigs->seqLen; i++){
            if(theAssembly->contigs->seq[i] != 'A' && theAssembly->contigs->seq[i] != 'T' && theAssembly->contigs->seq[i] != 'C' && theAssembly->contigs->seq[i] != 'G'){
                printf("Found an error in the assembly: theAssembly->contigs->seq[%i] = %c\n", i, theAssembly->contigs->seq[i]);
                return 0;
            }
        }
    }
    return 1;
}