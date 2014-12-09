#include "sources/neighbors.h"
#include <stdlib.h>
#include <math.h>
#include "sources/utils/algebra.h"
#include <cgds/BufferTop.h>
#include "sources/utils/boolean.h"

// evaluate distance between M[i,] and M[ii,]
double getDistance(double* M, int i, int ii, int ncol, double alpha, 
	bool simpleDists)
{
	// if simpleDists is ON, it means we are in stage 2 after convex optimization
	// ==> use full data, we know that now there are no NA's.
	if (simpleDists)
		return distance2(M+i*ncol, M+ii*ncol, ncol);

	// get distance for values per year
	double dist1 = 0.0;
	int valCount = 0; // number of not-NA fields
	int nobs = ncol-2; // ncol is 9+2 for our initial dataset (2001 to 2009)
	for (int j=0; j<nobs; j++)
	{
		if (!isnan(M[i*ncol+j]) && !isnan(M[ii*ncol+j]))
		{
			double diff = M[i*ncol+j] - M[ii*ncol+j];
			dist1 += diff*diff;
			valCount++;
		}
	}
	if (valCount > 0)
		dist1 /= valCount;

	// get distance for coordinates values
	double dist2 = 0.0;
	for (int j=nobs; j<ncol; j++)
	{
		double diff = M[i*ncol+j] - M[ii*ncol+j];
		dist2 += diff*diff;
	}
	dist2 /= 2.0; //to harmonize with above normalization
	if (valCount == 0)
		return sqrt(dist2); // no other choice

	//NOTE: adaptive alpha, the more NA's in vector, 
	//	  the more geo. coords. are taken into account
	alpha = (alpha >= 0.0 ? alpha : (double)valCount/nobs);
	return sqrt(alpha*dist1 + (1.0-alpha)*dist2);
}

// symmetrize neighborhoods lists (augmenting or reducing)
void symmetrizeNeighbors(List** neighborhoods, int nrow, int gmode)
{
	IndDist curNeighbI, curNeighbJ;
	for (int i=0; i<nrow; i++)
	{
		ListIterator* iterI = list_get_iterator(neighborhoods[i]);
		while (listI_has_data(iterI))
		{
			listI_get(iterI, curNeighbI);
			// check if neighborhoods[curNeighbI->index] has i
			bool reciproc = S_FALSE;
			List* neighbsJ = neighborhoods[curNeighbI.index];
			ListIterator* iterJ = list_get_iterator(neighbsJ);
			while (listI_has_data(iterJ))
			{
				listI_get(iterJ, curNeighbJ);
				if (curNeighbJ.index == i)
				{
					reciproc = S_TRUE;
					break;
				}
				listI_move_next(iterJ);
			}

			if (!reciproc)
			{
				if (gmode == 1)
				{
					// augmenting:
					// add (previously) non-mutual neighbor to neighbsJ
					list_insert_back(neighbsJ, i);
				}
				// test list_size() >= 2 because we don't allow empty neighborhoods
				else if (gmode == 0 && list_size(neighborhoods[i]) >= 2)
				{
					// reducing:
					// remove non-mutual neighbor to neighbsI
					listI_remove(iterI,BACKWARD);
				}
			}
			listI_move_next(iterI);
			listI_destroy(iterJ);
		}
		listI_destroy(iterI);
	}
}

// restrain neighborhoods: choose one per quadrant (for convex optimization)
void restrainToQuadrants(List** neighborhoods, int nrow, int ncol, double* M)
{
	IndDist curNeighbI;
	for (int i=0; i<nrow; i++)
	{
		ListIterator* iter = list_get_iterator(neighborhoods[i]);
		// choose one neighbor in each quadrant (if available); 
		// WARNING: multi-constraint optimization,
		//   > as close as possible to angle bissectrice
		//   > not too far from current data point

		// resp. SW,NW,SE,NE "best" neighbors :
		int bestIndexInDir[4] = {-1,-1,-1,-1};
		// corresponding "performances" :
		double bestPerfInDir[4] = {INFINITY,INFINITY,INFINITY,INFINITY};
		while (listI_has_data(iter))
		{
			listI_get(iter, curNeighbI);
			// get delta_x and delta_y to know in which quadrant 
			// we are and then get "index performance"
			// ASSUMPTION: all sites are distinct
			double deltaX = 
				M[i*ncol+(ncol-2)] - M[curNeighbI.index*ncol+(ncol-2)];
			double deltaY = 
				M[i*ncol+(ncol-1)] - M[curNeighbI.index*ncol+(ncol-1)];
			double angle = fabs(atan(deltaY/deltaX));
			// naive combination; [TODO: improve]
			double perf = curNeighbI.distance + fabs(angle-M_PI_4);
			// map {-1,-1} to 0, {-1,1} to 1 ...etc :
			int index = 2*(deltaX>0)+(deltaY>0);
			if (perf < bestPerfInDir[index])
			{
				bestIndexInDir[index] = curNeighbI.index;
				bestPerfInDir[index] = perf;
			}
			listI_move_next(iter);
		}

		// restrain neighborhood to the "best directions" found
		listI_reset_head(iter);
		while (listI_has_data(iter))
		{
			listI_get(iter, curNeighbI);
			// test list_size() <= 1 because we don't allow empty neighborhoods
			if (list_size(neighborhoods[i]) <= 1 ||
				curNeighbI.index==bestIndexInDir[0] || 
				curNeighbI.index==bestIndexInDir[1] || 
				curNeighbI.index==bestIndexInDir[2] || 
				curNeighbI.index==bestIndexInDir[3])
			{
				// OK, keep it
				listI_move_next(iter);
				continue;
			}
			// remove current node
			listI_remove(iter,FORWARD);
		}
		listI_destroy(iter);
	}
}

// Function to obtain neighborhoods.
// NOTE: alpha = weight parameter to compute distances; -1 means "adaptive"
List** getNeighbors_core(double* M, double alpha, int k, int gmode, 
	bool simpleDists, int nrow, int ncol)
{
	// prepare list buffers to get neighborhoods
	// (OK for small to moderate values of k)
	BufferTop** bufferNeighbs = 
		(BufferTop**)malloc(nrow*sizeof(BufferTop*));
	for (int i=0; i<nrow; i++)
		bufferNeighbs[i] = buffertop_new(IndDist, k, MIN_T, 2);

	// MAIN LOOP

	// for each row in M, find its k nearest neighbors
	for (int i=0; i<nrow; i++)
	{
		// for each potential neighbor...
		for (int ii=0; ii<nrow; ii++)
		{
			if (ii == i) 
				continue;

			// evaluate distance from M[i,] to M[ii,]
			double distance = 
				getDistance(M, i, ii, ncol, alpha, simpleDists);

			// (try to) add index 'ii' + distance to bufferNeighbs[i]
			IndDist id = (IndDist){.index=ii, .distance=distance};
			buffertop_tryadd(bufferNeighbs[i], id, distance);
		}
	}

	// free buffers and transfer their contents into lists easier to process
	List** neighborhoods = (List**)malloc(nrow*sizeof(List*));
	for (int i=0; i<nrow; i++)
	{
		neighborhoods[i] = buffertop_2list(bufferNeighbs[i]);
		buffertop_destroy(bufferNeighbs[i]);
	}
	free(bufferNeighbs);

	// OPTIONAL MUTUAL KNN
	if (gmode==0 || gmode==1)
	{
		// additional processing to symmetrize neighborhoods (augment or not)
		symmetrizeNeighbors(neighborhoods, nrow, gmode);
	}
	else if (gmode==3)
	{
		// choose one neighbor per quadrant (for convex optimization)
		restrainToQuadrants(neighborhoods, nrow, ncol, M);
	}
	// nothing to do if gmode==2 (simple assymmetric kNN)

	return neighborhoods;
}
