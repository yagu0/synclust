#include "tests/helpers.h"
#include "sources/dijkstra.h"

//custom graph 1 (connected!)
void test_dijkstra1()
{
	int n = 10;
	double distsIn[100] = 
	{
		NAN,4.0,6.0,3.0,NAN,1.0,NAN,NAN,NAN,NAN,
		4.0,NAN,NAN,5.0,NAN,1.0,3.0,NAN,1.0,NAN,
		6.0,NAN,NAN,NAN,NAN,7.0,1.0,NAN,NAN,NAN,
		3.0,5.0,NAN,NAN,4.0,NAN,NAN,NAN,NAN,1.0,
		NAN,NAN,NAN,4.0,NAN,NAN,NAN,2.0,3.0,NAN,
		1.0,1.0,7.0,NAN,NAN,NAN,NAN,3.0,NAN,NAN,
		NAN,3.0,1.0,NAN,NAN,NAN,NAN,NAN,1.0,NAN,
		NAN,NAN,NAN,NAN,2.0,3.0,NAN,NAN,NAN,1.0,
		NAN,1.0,NAN,NAN,3.0,NAN,1.0,NAN,NAN,NAN,
		NAN,NAN,NAN,1.0,NAN,NAN,NAN,1.0,NAN,NAN
	};
	double* distsOut;

	distsOut = dijkstra_core(distsIn, 0, 10);
	//as by-hand computed, distances should be as follow
	double shouldOutput0[10] = {0.0,2.0,5.0,3.0,6.0,1.0,4.0,4.0,3.0,4.0};
	ASSERT_TRUE(checkEqualV(shouldOutput0, distsOut, n));
	free(distsOut);

	distsOut = dijkstra_core(distsIn, 7, 10);
	//as by-hand computed, distances should be as follow
	double shouldOutput7[10] = {4.0,4.0,7.0,2.0,2.0,3.0,6.0,0.0,5.0,1.0};
	ASSERT_TRUE(checkEqualV(shouldOutput7, distsOut, n));
	free(distsOut);
}

//custom graph 2 (connected!)
void test_dijkstra2()
{
	int n = 10;
	// same as graph 1 above, but link between 1 and 5 is now 4.0 instead of 1.0
	double distsIn[100] = 
	{
		NAN,4.0,6.0,3.0,NAN,1.0,NAN,NAN,NAN,NAN,
		4.0,NAN,NAN,5.0,NAN,4.0,3.0,NAN,1.0,NAN,
		6.0,NAN,NAN,NAN,NAN,7.0,1.0,NAN,NAN,NAN,
		3.0,5.0,NAN,NAN,4.0,NAN,NAN,NAN,NAN,1.0,
		NAN,NAN,NAN,4.0,NAN,NAN,NAN,2.0,3.0,NAN,
		1.0,4.0,7.0,NAN,NAN,NAN,NAN,3.0,NAN,NAN,
		NAN,3.0,1.0,NAN,NAN,NAN,NAN,NAN,1.0,NAN,
		NAN,NAN,NAN,NAN,2.0,3.0,NAN,NAN,NAN,1.0,
		NAN,1.0,NAN,NAN,3.0,NAN,1.0,NAN,NAN,NAN,
		NAN,NAN,NAN,1.0,NAN,NAN,NAN,1.0,NAN,NAN
	};
	double* distsOut;

	distsOut = dijkstra_core(distsIn, 0, 10);
	//as by-hand computed, distances should be as follow
	double shouldOutput0[10] = {0.0,4.0,6.0,3.0,6.0,1.0,6.0,4.0,5.0,4.0};
	ASSERT_TRUE(checkEqualV(shouldOutput0, distsOut, n));
	free(distsOut);

	distsOut = dijkstra_core(distsIn, 7, 10);
	//as by-hand computed, distances should be as follow
	double shouldOutput7[10] = {4.0,6.0,7.0,2.0,2.0,3.0,6.0,0.0,5.0,1.0};
	ASSERT_TRUE(checkEqualV(shouldOutput7, distsOut, n));
	free(distsOut);
}
