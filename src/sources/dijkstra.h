#ifndef SYNCLUST_DIJKSTRA_H
#define SYNCLUST_DIJKSTRA_H

// Dijkstra from index start : return vector of distances to every other vertex
double* dijkstra_core(
	double* pDistsIn, 
	int start, 
	int n
);

#endif
