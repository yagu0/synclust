#include "tests/helpers.h"
#include "sources/connexity.h"

//completely disconnected graph (no edges)
void test_connexity1()
{
	int n = 10;
	int* lengthNIix = (int*)calloc(n,sizeof(int));
	int** NIix = (int**)malloc(n*sizeof(int*));
	for (int i=0; i<n; i++) NIix[i] = NULL;
	int* cc = getConnectedComponents_core(NIix, lengthNIix, n);
	//cc should contain all integers from 1 to n
	ASSERT_TRUE(countDistinctValues(cc,n) == n);
	free(lengthNIix);
	free(NIix);
	free(cc);
}

//bipartite graph
void test_connexity2()
{
	int n = 10;
	int* lengthNIix = (int*)malloc(n*sizeof(int));
	int** NIix = (int**)malloc(n*sizeof(int*));
	for (int i=0; i<n; i++)
	{
		lengthNIix[i] = 1;
		NIix[i] = (int*)malloc(sizeof(int));
	}
	// connect 0 with 1, 2 with 3 ...
	for (int i=0; i<n; i+=2)
	{
		NIix[i][0] = i+1;
		NIix[i+1][0] = i;
	}
	int* cc = getConnectedComponents_core(NIix, lengthNIix, n);
	//cc should contain exactly n/2 integers
	ASSERT_TRUE(countDistinctValues(cc,n) == n/2);
	free(lengthNIix);
	for (int i=0; i<n; i++)
		free(NIix[i]);
	free(NIix);
	free(cc);
}

//~ //cyclic graph
//~ void test_connexity3() {
	//~ int n = 10;
	//~ int* adjMat = (int*)calloc(n*n,sizeof(int));
	//~ // connect 0 with 1, 1 with 2 ...
	//~ for (int i=0; i<n; i++) {
		//~ adjMat[i+n*((i+1)%n)] = TRUE;
		//~ adjMat[(i+1)%n+n*i] = TRUE;
	//~ }
	//~ int* cc = getConnectedComponents_core(adjMat, n);
	//~ //cc should contain only one value (presumably '1')
	//~ warnIfFails(countDistinctValues(cc,n) == 1, "c3", "");
	//~ free(adjMat);
	//~ free(cc);
//~ }
//~ 
//~ //custom graph with 3 connex components
//~ void test_connexity4() {
	//~ int n = 10;
	//~ 
	//~ 
		//~ {2,4},
		//~ {2,4},
		//~ {0,1},
		//~ {5,8,9},
		//~ {0,1},
		//~ {3,7},
		//~ {},
		//~ {5,9},
		//~ {3},
		//~ {3,7}
	//~ 
	//~ int adjMat[100] = 
	//~ {
		//~ 0,0,1,0,1,0,0,0,0,0,
		//~ 0,0,1,0,1,0,0,0,0,0,
		//~ 1,1,0,0,0,0,0,0,0,0,
		//~ 0,0,0,0,0,1,0,0,1,1,
		//~ 1,1,0,0,0,0,0,0,0,0,
		//~ 0,0,0,1,0,0,0,1,0,0,
		//~ 0,0,0,0,0,0,0,0,0,0,
		//~ 0,0,0,0,0,1,0,0,0,1,
		//~ 0,0,0,1,0,0,0,0,0,0,
		//~ 0,0,0,1,0,0,0,1,0,0
	//~ };
	//~ int* cc = getConnectedComponents_core(adjMat, n);
	//~ //cc should contain exactly 3 values
	//~ warnIfFails(countDistinctValues(cc,n) == 3, "c4", "");
	//~ free(cc);
//~ }
//~ 
//~ //custom graph, 1 connex component
//~ void test_connexity5() {
	//~ int n = 10;
	//~ int adjMat[100] = 
	//~ {
		//~ 0,0,1,1,0,0,0,1,0,0,
		//~ 0,0,1,0,1,0,1,0,0,0,
		//~ 1,1,0,0,0,0,0,0,0,0,
		//~ 1,0,0,0,0,1,0,0,1,1,
		//~ 0,1,0,0,0,0,0,0,0,0,
		//~ 0,0,0,1,0,0,0,1,0,0,
		//~ 0,1,0,0,0,0,0,0,0,0,
		//~ 1,0,0,0,0,1,0,0,0,1,
		//~ 0,0,0,1,0,0,0,0,0,0,
		//~ 0,0,0,1,0,0,0,1,0,0
	//~ };
	//~ int* cc = getConnectedComponents_core(adjMat, n);
	//~ //cc should contain only one value (presumably '1')
	//~ warnIfFails(countDistinctValues(cc,n) == 1, "c5", "");
	//~ free(cc);
//~ }
