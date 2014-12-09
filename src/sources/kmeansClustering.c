#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "sources/utils/boolean.h"
#include "sources/kmeansClustering.h"

// auxiliary function to obtain a random sample of 1..n with K elements
void sample(int* centers, int n, int K)
{
	// refVect = (0,1,...,n-1,n)
	int* refVect = (int*)malloc(n*sizeof(int));
	for (int i=0; i<n; i++)
		refVect[i] = i;

	int curSize = n; // current size of the sampling set
	for (int j=0; j<K; j++)
	{
		// pick an index in sampling set:
		int index = rand()%curSize;
		centers[j] = refVect[index];
		// move this index outside of sampling set:
		refVect[index] = refVect[--curSize];
	}

	free(refVect);
}

// auxiliary function to compare two sets of centers
int unequalCenters(int* ctrs1, int* ctrs2, int n, int K)
{
	// HACK: special case at initialization, ctrs2 = 0
	if (K > 1 && ctrs2[0]==0 && ctrs2[1]==0)
		return S_TRUE;

	// compVect[i] == 1 iff index i is found in ctrs1 or ctrs2
	int* compVect = (int*)calloc(n,sizeof(int));

	int kountNonZero = 0; // count non-zero elements in compVect
	for (int j=0; j<K; j++)
	{
		if (compVect[ctrs1[j]] == 0)
			kountNonZero++;
		compVect[ctrs1[j]] = 1;
		if (compVect[ctrs2[j]] == 0)
			kountNonZero++;
		compVect[ctrs2[j]] = 1;
	}

	free(compVect);

	// if we found more than K non-zero elements, ctrs1 and ctrs2 differ
	return (kountNonZero > K);
}

// assign a vector (represented by its distances to others, as distances[index,])
// to a cluster, represented by its center ==> output is integer in 0..K-1
int assignCluster(int index, double* distances, int* centers, int n, int K)
{
	int minIndex = 0;
	double minDist = distances[index*n+centers[0]];

	for (int j=1; j<K; j++)
	{
		if (distances[index*n+centers[j]] < minDist)
		{
			minDist = distances[index*n+centers[j]];
			minIndex = j;
		}
	}

	return minIndex;
}

// k-means based on a distance matrix (nstart=10, maxiter=100)
int* kmeansWithDistances_core(
	double* distances, int n, int K, int nstart, int maxiter)
{
	int* bestClusts = (int*)malloc(n*sizeof(int));
	double bestDistor = INFINITY;
	double avgClustSize = (double)n/K;
	int* ctrs = (int*)malloc(K*sizeof(int));
	int* oldCtrs = (int*)malloc(K*sizeof(int));
	Vector** clusters = (Vector**)malloc(K*sizeof(Vector*));
	for (int j=0; j<K; j++)
		clusters[j] = vector_new(int);

	// set random number generator seed
	srand(time(NULL));

	for (int startKount=0; startKount < nstart; startKount++)
	{
		// centers (random) [re]initialization
		sample(ctrs,n,K);
		for (int j=0; j<K; j++)
			oldCtrs[j] = 0;
		int kounter = 0;

		/*************
		 *  main loop
		 *************/

		// while old and "new" centers differ..
		while (unequalCenters(ctrs,oldCtrs,n,K) && kounter++ < maxiter)
		{
			// (re)initialize clusters to empty sets
			for (int j=0; j<K; j++)
				vector_clear(clusters[j]);

			// estimate clusters belongings
			for (int i=0; i<n; i++)
			{
				int affectation = assignCluster(i, distances, ctrs, n, K);
				vector_push(clusters[affectation], i);
			}

			// copy current centers to old centers
			for (int j=0; j<K; j++)
				oldCtrs[j] = ctrs[j];

			// recompute centers
			for (int j=0; j<K; j++)
			{
				int minIndex = -1;
				double minSumDist = INFINITY;
				VectorIterator* iter1 = vector_get_iterator(clusters[j]);
				vectorI_reset_begin(iter1);
				while (vectorI_has_data(iter1))
				{
					int index1; vectorI_get(iter1, index1);
					// attempt to use current index as center
					double sumDist = 0.0;
					VectorIterator* iter2 = vector_get_iterator(clusters[j]);
					vectorI_reset_begin(iter2);
					while (vectorI_has_data(iter2))
					{
						int index2; vectorI_get(iter2, index2);
						sumDist += distances[index1*n+index2];
						vectorI_move_next(iter2);
					}
					if (sumDist < minSumDist)
					{
						minSumDist = sumDist;
						minIndex = index1;
					}
					vectorI_destroy(iter2);
					vectorI_move_next(iter1);
				}
				if (minIndex >= 0)
					ctrs[j] = minIndex;
				// HACK: some 'random' index (a cluster should never be empty)
				// this case should never happen anyway 
				// (pathological dataset with replicates)
				else
					ctrs[j] = 0;
				vectorI_destroy(iter1);
			}
		} /***** end main loop *****/

		// finally compute distorsions :
		double distor = 0.0;
		for (int j=0; j<K; j++)
		{
			VectorIterator* iter = vector_get_iterator(clusters[j]);
			vectorI_reset_begin(iter);
			while (vectorI_has_data(iter))
			{
				int index; vectorI_get(iter, index);
				distor += distances[index*n+ctrs[j]];
				vectorI_move_next(iter);
			}
			vectorI_destroy(iter);
		}
		if (distor < bestDistor)
		{
			// copy current clusters into bestClusts
			for (int j=0; j<K; j++)
			{
				VectorIterator* iter = vector_get_iterator(clusters[j]);
				vectorI_reset_begin(iter);
				while (vectorI_has_data(iter))
				{
					int index; vectorI_get(iter, index);
					bestClusts[index] = j;
					vectorI_move_next(iter);
				}
				vectorI_destroy(iter);
			}
			bestDistor = distor;
		}
	}

	free(ctrs);
	free(oldCtrs);
	for (int j=0; j<K; j++)
		vector_destroy(clusters[j]);
	free(clusters);

	return bestClusts;
}
