#include "sources/dijkstra.h"
#include <stdlib.h>
#include "sources/utils/boolean.h"
#include <math.h>

// Dijkstra from index start : return vector of distances to every other vertex
// TODO: use a good priority queue, and pass NI instead of pDistsIn (+ linear preprocessing)
double* dijkstra_core(double* pDistsIn, int start, int n) {

	// initalisations
	double* pDistsOut = (double*)malloc(n*sizeof(double));
	bool* visited = (bool*)malloc(n*sizeof(bool));
	for (int i=0; i<n; i++) {
		pDistsOut[i] = INFINITY;
		visited[i] = S_FALSE; // nothing seen so far
	}
	pDistsOut[start] = 0.0; // we are at distance 0 from self

	while (S_TRUE)
	{
		double minGeodDist = INFINITY;

		// n1 <-- node in "unseen" with smallest geodesic distance
		// NOTE: on first loop, n1 == start
		int n1 = 0;
		for (int i=0; i<n; i++)
		{
			if (!visited[i] && pDistsOut[i] < minGeodDist)
			{
				n1 = i;
				minGeodDist = pDistsOut[i];
			}
		}

		if (minGeodDist == INFINITY)
			break; // all is explored

		visited[n1] = S_TRUE; // n1 is expanded

		// For n2 running through neighbors of n1
		for (int n2 = 0; n2<n; n2++)
		{
			int ind_n12 = n1*n+n2; // or n1+n*n2 (symmetry)
			if (!isnan(pDistsIn[ind_n12]))
			{
				// check if we'd better go through n1 (on the way from start to n2)
				if (pDistsOut[n2] > pDistsOut[n1] + pDistsIn[ind_n12])
					pDistsOut[n2] = pDistsOut[n1] + pDistsIn[ind_n12];
			}
		}
	}
	free(visited);

	return pDistsOut;
}
