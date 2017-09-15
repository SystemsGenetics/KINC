#ifndef _MEAN_SHIFT_
#define _MEAN_SHIFT_


#include "PairWiseCluster.h"

// Function to peform the clustering
PairWiseClusters * clustering(float *a2, int x, float *b2, int y, int n2,
    EMatrix * ematrix, CCMParameters params, float bw, int level);


#endif
