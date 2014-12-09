#include "tests/helpers.h"
#include "sources/kmeansClustering.h"
#include <math.h>

//easy data: already clustered, one cluster = 1 vertex in equilateral triangle
void test_kmeansClustering1()
{
	int n=99;
	double* distances = (double*)malloc(n*n*sizeof(double));
	for (int i=0; i<n*n; i++)
		distances[i] = 1.0;
	int clustCount = 3, clustSize = n/clustCount; //33

	for (int kounter=0; kounter<clustCount; kounter++)
	{
		//cluster k: kounter*33...(kounter+1)*33-1
		for (int i=kounter*clustSize; i<(kounter+1)*clustSize; i++)
		{
			for (int j=kounter*clustSize; j<(kounter+1)*clustSize; j++)
				distances[i+n*j] = 0.0; //high-density cluster...
		}
	}

	//call to clustering algorithm
	int* clusters = kmeansWithDistances_core(distances, n, clustCount, 10, 100);

	ASSERT_TRUE(countDistinctValues(clusters, n) == clustCount);
	ASSERT_TRUE(checkClustersProportions(clusters, n, clustCount, 1e-10));
	free(distances);
	free(clusters);
}

//three isotropic (well separated) gaussian clusters
void test_kmeansClustering2()
{
	// generate 2D data
	int n=99, d=2;
	double* M = (double*)malloc(n*d*sizeof(double));
	int clustCount = 3, clustSize = n/clustCount; //33

	double ctrs[3][2] =
	{
		{-3.0,-3.0},
		{0.0,0.0},
		{3.0,3.0}
	};

	srand(time(NULL));
	for (int kounter=0; kounter<clustCount; kounter++)
	{
		//cluster k: kounter*33...(kounter+1)*33-1
		for (int i=kounter*clustSize; i<(kounter+1)*clustSize; i++)
		{
			double U = (double)rand()/RAND_MAX;
			double V = (double)rand()/RAND_MAX;
			double fact = sqrt(-2*log(U));
			M[i+n*0] = ctrs[kounter][0] + fact * cos(2*M_PI*V);
			M[i+n*1] = ctrs[kounter][1] + fact * sin(2*M_PI*V);
		}
	}

	// compute distances matrix
	double* distances = (double*)calloc(n*n,sizeof(double));
	for (int i=0; i<n; i++)
	{
		for (int ii=0; ii<n; ii++)
		{
			double distance = 0.0;
			for (int j=0; j<d; j++)
				distance += (M[i+n*j] - M[ii+n*j])*(M[i+n*j] - M[ii+n*j]);
			distances[i+n*ii] = sqrt(distance);
		}
	}
	free(M); //no need for initial data anymore

	//call to clustering algorithm
	int* clusters = kmeansWithDistances_core(distances, n, clustCount, 10, 100);

	ASSERT_TRUE(countDistinctValues(clusters, n) == clustCount);
	//test that each cluster accounts for 1/3 of total data, +/- 10%
	ASSERT_TRUE(checkClustersProportions(clusters, n, clustCount, 0.1));
	free(distances);
	free(clusters);
}
