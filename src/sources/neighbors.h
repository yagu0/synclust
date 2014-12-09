#ifndef SYNCLUST_NEIGHBORS_H
#define SYNCLUST_NEIGHBORS_H

#include "sources/utils/boolean.h"
#include <cgds/List.h>

// evaluate distance between M[i,] and M[ii,]
double getDistance(
	double* M, 
	int i, 
	int ii, 
	int ncol, 
	double alpha, 
	bool simpleDists
);

// symmetrize neighborhoods lists (augmenting or reducing)
void symmetrizeNeighbors(
	List** neighborhoods, 
	int nrow, 
	int gmode
);

// restrain neighborhoods: choose one per quadrant (for convex optimization)
void restrainToQuadrants(
	List** neighborhoods, 
	int nrow, 
	int ncol, 
	double* M
);

// structure to store a neighbor index and the distance to this neighbor
typedef struct IndDist {
	int index;
	double distance;
} IndDist;

// Function to obtain neighborhoods.
// NOTE: alpha = weight parameter to compute distances; -1 means "adaptive"
List** getNeighbors_core(
	double* M, 
	double alpha, 
	int k, 
	int gmode, 
	bool simpleDists, 
	int nrow, 
	int ncol
);

#endif
