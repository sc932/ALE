#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const static double lnfactconst2 = 0.918938533204672741780329;
const static double minLogLike = -120.0;

//uses Stirlings approximation to high precision
double lnfact2(double input){
  return (input - 0.5)*log(input) - input + lnfactconst2 - 1.0/(12.0*input) - 1.0/(360.0*input*input*input) - 1.0/(1260.0*input*input*input*input*input);
}

// finds the poisson pmf at value k for mean lambda
double poissonPMF(double k, double lambda){
    return k*log(lambda) - lambda - lnfact2(k + 1);
}

int main(int argc, char **argv){
    if (argc < 4) {
        printf("Usage: %s ALEoutput1 ALEoutput2 newALEoutput\n(C) 2011 Scott Clark\nTakes two ALE output files and combines them into a new output\n\n", argv[0]);
        return 0;
    }
    
    // attempt to open the first input file
    FILE *ins1 = fopen(argv[argc - 3], "r");
    if(ins1 == NULL){
        printf("Error! Could not open the first input file: %s\n", argv[argc - 3]);
    }
    
    // attempt to open the second input file
    FILE *ins2 = fopen(argv[argc - 2], "r");
    if(ins2 == NULL){
        printf("Error! Could not open the second input file: %s\n", argv[argc - 2]);
    }
    
    // attempt to open the output file
    FILE *fout = fopen(argv[argc - 1], "w");
    if(fout == NULL){
        printf("Error! Could not open the output file: %s\n", argv[argc - 1]);
    }
    
    int keepGoing = 1, i;
    char name1[256], name2[256];
    char temp0[10], temp1[256], temp2[256], temp3[256], temp4[256], temp5[256], temp6[256];
    int len1, len2;
    double depth1, depth2;
    double placement1, placement2;
    double kmer;
    double depthProb, totalProb, placeProb;
    double depthTotal;
    while(keepGoing > 0){
      // read in the header
      keepGoing = fscanf( ins1, "%255s%i%10s%255s%255s%255s%255s%255s", name1, &len1, temp1, temp2, temp3, temp4, temp5, temp6);
      keepGoing = fscanf( ins2, "%255s%i%10s%255s%255s%255s%255s%255s", name2, &len2, temp1, temp2, temp3, temp4, temp5, temp6);
      //printf("%s %i %s %s %s %s %s\n", name1, len1, temp1, temp2, temp3, temp4, temp5, temp6);
      //printf("%s %i %s %s %s %s %s\n", name2, len2, temp1, temp2, temp3, temp4, temp5, temp6);
      if(strcmp(name1, name2) != 0 || len1 != len2){
	  printf("Error! The contigs in the ALE outputs do not match! Make sure they are generated using the same assembly!\n");
	  printf("Found name1: %s, name2: %s, len1: %i, len2: %i.\n", name1, name2, len1, len2);
	  return 0;
      }
      if(keepGoing == 0){
	  break;
      }
      for(i = 0; i < len1; i++){
	keepGoing = fscanf( ins1, "%lf%lf%lf%lf%lf", &depth1, &depthProb, &placement1, &kmer, &totalProb);
	keepGoing = fscanf( ins2, "%lf%lf%lf%lf%lf", &depth2, &depthProb, &placement2, &kmer, &totalProb);
	//printf("1: %lf %lf %lf %lf %lf\n", depth1, depthProb, placement1, kmer, totalProb);
	//printf("2: %lf %lf %lf %lf %lf\n", depth2, depthProb, placement2, kmer, totalProb);
	depthTotal += depth1 + depth2;
	//printf("Depth avg (%i) = %lf\n", len1, depthTotal);
      }
    }
    
    depthTotal = depthTotal/((float)len1);
    printf("Depth avg (%i) = %lf\n", len1, depthTotal);
    
    fclose(ins1);
    fclose(ins2);
    
    // attempt to reopen the first input file
    ins1 = fopen(argv[argc - 3], "r");
    if(ins1 == NULL){
        printf("Error! Could not reopen the first input file: %s\n", argv[argc - 3]);
    }
    
    // attempt to reopen the second input file
    ins2 = fopen(argv[argc - 2], "r");
    if(ins2 == NULL){
        printf("Error! Could not reopen the second input file: %s\n", argv[argc - 2]);
    }
    
    keepGoing = 1;
    while(keepGoing > 0){
      // read in the header
      keepGoing = fscanf( ins1, "%255s%i%10s%255s%255s%255s%255s%255s", name1, &len1, temp1, temp2, temp3, temp4, temp5, temp6);
      keepGoing = fscanf( ins2, "%255s%i%10s%255s%255s%255s%255s%255s", name2, &len2, temp1, temp2, temp3, temp4, temp5, temp6);
      if(strcmp(name1, name2) != 0 || len1 != len2){
	  printf("Error! The contigs in the ALE outputs do not match! Make sure they are generated using the same assembly!\n");
	  printf("Found name1: %s, name2: %s, len1: %i, len2: %i.\n", name1, name2, len1, len2);
	  return 0;
      }
      //printf("Keep going: %i\n", keepGoing);
      if(keepGoing < 1){
	  break;
      }else{
	fprintf(fout, "%s %i %s %s %s %s %s %s\n", name1, len1, temp1, temp2, temp3, temp4, temp5, temp6);
	for(i = 0; i < len1; i++){
	  keepGoing = fscanf( ins1, "%lf%lf%lf%lf%lf", &depth1, &depthProb, &placement1, &kmer, &totalProb);
	  keepGoing = fscanf( ins2, "%lf%lf%lf%lf%lf", &depth2, &depthProb, &placement2, &kmer, &totalProb);
	  depthProb = poissonPMF(depth1 + depth2, depthTotal);
	  if(depthProb < minLogLike){
	      depthProb = minLogLike;
	  }
	  if(depth1 != 0.0 && depth2 != 0.0){
	    placeProb = log((depth1*exp(placement1) + depth2*exp(placement2))/(depth1 + depth2));
	  }else if(depth1 != 0.0 && depth2 == 0.0){
	    placeProb = placement1;
	  }else if(depth2 != 0.0 && depth1 == 0.0){
	    placeProb = placement2;
	  }else{
	    placeProb = minLogLike;
	  }
	  //printf("%lf %lf %lf %lf %lf\n", depth1 + depth2, depthProb, placeProb, kmer, depthProb+ placeProb+ kmer);
	  fprintf(fout, "%lf %lf %lf %lf %lf\n", depth1 + depth2, depthProb, placeProb, kmer, depthProb+ placeProb+ kmer);
	}
      }
    }
    fclose(ins1);
    fclose(ins2);
    fclose(fout);

//     while(keepGoing > 0){
//         // the is the contig header
//         keepGoing = fscanf( ins1, "%255s%i%10s%255s%255s%255s%255s%255s", name1, &len1, temp0, temp1, temp2, temp3, temp4, temp5);
//         keepGoing = fscanf( ins2, "%255s%i%10s%255s%255s%255s%255s%255s", name2, &len2, temp0, temp1, temp2, temp3, temp4, temp5);
//         if(strcmp(name1, name2) != 0 || len1 != len2){
//             printf("Error! The contigs in the ALE outputs do not match! Make sure they are generated using the same assembly!\n");
//             return 0;
//         }
//         if(keepGoing == 0){
//             break;
//         }
//         printf("%s %i %s %s %s %s %s %s\n", name2, len2, temp0, temp1, temp2, temp3, temp4, temp5);
//         depthTotal = 0.0;
//         depth1 = malloc(len1*sizeof(double));
//         depth2 = malloc(len1*sizeof(double));
//         placement1 = malloc(len1*sizeof(double));
//         placement2 = malloc(len1*sizeof(double));
//         kmer = malloc(len1*sizeof(double));
//         for(i = 0; i < len1; i++){
//             keepGoing = fscanf( ins1, "%lf%lf%lf%lf%lf", &depth1[i], &depthProb, &placement1[i], &kmer[i], &totalProb);
//             keepGoing = fscanf( ins2, "%lf%lf%lf%lf%lf", &depth2[i], &depthProb, &placement2[i], &kmer[i], &totalProb);
//             //printf("in1: %lf %lf %lf %lf %lf\n", depth1[i], depthProb, placement1[i], kmer[i], totalProb);
//             //printf("in2: %lf %lf %lf %lf %lf\n", depth2[i], depthProb, placement2[i], kmer[i], totalProb);
//             depthTotal += depth1[i] + depth2[i];
//         }
//         depthTotal = depthTotal/(float)len1;
//         fprintf(fout, ">%s %i > depth : ln(depthLike) : ln(placeLike) : ln(kmerLike) : ln(totalLike)\n", name1, len1);
//         for(i = 0; i < len1; i++){
//             depthProb = poissonPMF(depth1[i] + depth2[i], depthTotal);
//             if(depthProb < minLogLike){
//                 depthProb = minLogLike;
//             }
//             fprintf(fout, "%lf %lf %lf %lf %lf\n", depth1[i] + depth2[i], depthProb, log((depth1[i]*exp(placement1[i]) + depth2[i]*exp(placement2[i]))/(depth1[i] + depth2[i])), kmer[i], depthProb + log((depth1[i]*exp(placement1[i]) + depth2[i]*exp(placement2[i]))/(depth1[i] + depth2[i])) + kmer[i]);
//         }
//     }
    
}