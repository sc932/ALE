// geneTree.h

#ifndef _GENE_TREE_H_
#define _GENE_TREE_H_

#include <stdlib.h>
#include <stdio.h>
#include "ALEhelpers.h"

#define START_LENGTH        5

struct LeafNode{
  int length;
  int current;
  int indicies[][2];
};

struct TreeBranch{
  struct TreeBranch *subBranches[4];// = {NULL,NULL,NULL,NULL};
  struct LeafNode *leaf;
};

typedef struct LeafNode leafNode_t;
typedef struct TreeBranch treeBranch_t;

int OutputIndicies(treeBranch_t *pRoot, const char sequence[], const int klen, int index[][2]);

treeBranch_t MakeTree(const assembly_t theAssem, const int klen, const int NUM_ASSEMBLY_PARTS);

int AddSeqToTree(const char sequence[], const int offset, const int klen, treeBranch_t *pRoot, const int assemPart);

#endif
