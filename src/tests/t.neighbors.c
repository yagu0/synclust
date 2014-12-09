#include "tests/helpers.h"
#include "sources/neighbors.h"
#include <cgds/List.h>

int** list2int(List** L, int n)
{
	int** I = (int**)malloc(n*sizeof(int*));
	for (int i=0; i<n; i++)
	{
		int listSize = list_size(L[i]);
		I[i] = (int*)malloc((1+listSize)*sizeof(int));
		I[i][0] = listSize;
		ListIterator* iterJ = list_get_iterator(L[i]);
		for (int j=1; j<=listSize; j++)
		{
			listI_get(iterJ, I[i][j]);
			listI_move_next(iterJ);
		}
		listI_destroy(iterJ);
	}
	return I;
}

//10 lines 12 columns, only NA except on antidiagonal
// ==> distances computed with coordinates only
void test_neighbors1()
{
	int n = 10, m=12;
	double M[120] = 
	{
		NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,1.0,0.0,0.0,
		NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,1.0,NAN,1.0,1.0,
		NAN,NAN,NAN,NAN,NAN,NAN,NAN,1.0,NAN,NAN,2.0,2.0,
		NAN,NAN,NAN,NAN,NAN,NAN,1.0,NAN,NAN,NAN,3.0,3.0,
		NAN,NAN,NAN,NAN,NAN,1.0,NAN,NAN,NAN,NAN,4.0,4.0,
		NAN,NAN,NAN,NAN,1.0,NAN,NAN,NAN,NAN,NAN,5.0,5.0,
		NAN,NAN,NAN,1.0,NAN,NAN,NAN,NAN,NAN,NAN,6.0,6.0,
		NAN,NAN,1.0,NAN,NAN,NAN,NAN,NAN,NAN,NAN,7.0,7.0,
		NAN,1.0,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,8.0,8.0,
		1.0,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,9.0,9.0
	};

	double alphas[4] = {-1.0, 0.0, 0.5, 1.0};
	int k = 2; // no need for more
	for (int j=0; j<4; j++)
	{
		double alpha = alphas[j]; // no impact
		for (int gmode=0; gmode<4; gmode++)
		{
			List** L = getNeighbors_core(M, alpha, k, gmode, S_FALSE, n, m);
			int** NIix = list2int(L, n);
			for (int jj=0; jj<n; jj++)
				list_destroy(L[jj]);
			free(L);
			for (int jj=1; jj<n-1; jj++)
			{
				// interior points: 2 neighbors, below and above [except for gmode==1 and jj==2 or n-2]
				if (gmode==1 && (jj==2 || jj==n-3))
				{
					ASSERT_TRUE(
						NIix[jj][0] == 3 && (NIix[jj][1] == jj-1 || NIix[jj][2] == jj-1 || NIix[jj][3] == jj-1));
					ASSERT_TRUE(
						NIix[jj][0] == 3 && (NIix[jj][1] == jj+1 || NIix[jj][2] == jj+1 || NIix[jj][3] == jj+1));
					if (jj==2)
					{
						ASSERT_TRUE(
							NIix[jj][0] == 3 && (NIix[jj][1] == jj-2 || NIix[jj][2] == jj-2 || NIix[jj][3] == jj-2));
					}
					else if (jj==n-3)
					{
						ASSERT_TRUE(
							NIix[jj][0] == 3 && (NIix[jj][1] == jj+2 || NIix[jj][2] == jj+2 || NIix[jj][3] == jj+2));
					}
				}
				else
				{
					// one neighb below, one neighb above
					ASSERT_TRUE(
						NIix[jj][0] == 2 && (NIix[jj][1] == jj-1 || NIix[jj][2] == jj-1));
					ASSERT_TRUE(
						NIix[jj][0] == 2 && (NIix[jj][1] == jj+1 || NIix[jj][2] == jj+1));
				}
			}
			// boundary points in mode 1 (augmented kNN) or 2 (normal kNN) also have 2 neighbors
			if (gmode==1 || gmode==2)
			{
				ASSERT_TRUE(NIix[0][0] == 2 && (NIix[0][1] == 1 || NIix[0][2] == 1));
				ASSERT_TRUE(NIix[0][0] == 2 && (NIix[0][1] == 2 || NIix[0][2] == 2));
				ASSERT_TRUE(NIix[n-1][0] == 2 && (NIix[n-1][1] == n-2 || NIix[n-1][2] == n-2));
				ASSERT_TRUE(NIix[n-1][0] == 2 && (NIix[n-1][1] == n-3 || NIix[n-1][2] == n-3));
			}
			else
			{
				// in mutual kNN or in 'by-quadrants' modes, they only have 1 neighbor.
				ASSERT_TRUE(NIix[0][0] == 1 && NIix[0][1] == 1);
				ASSERT_TRUE(NIix[n-1][0] == 1 && NIix[n-1][1] == n-2);
			}

			for (int i=0; i<n; i++)
				free(NIix[i]);
			free(NIix);
		}
	}
}

//10 lines, 6 columns, some NA's...
void test_neighbors2()
{
	int n = 6, m=10;
	double M[60] = 
	{
		NAN,1.0,0.0,NAN,0.0,NAN,3.0,NAN,0.0,0.0,
		1.0,2.0,0.0,NAN,NAN,0.0,2.0,3.0,1.0,1.0,
		2.0,3.0,0.0,NAN,NAN,1.0,1.0,2.0,2.0,2.0,
		3.0,NAN,0.0,NAN,1.0,NAN,NAN,1.0,3.0,3.0,
		2.0,3.0,NAN,3.0,1.0,2.0,1.0,NAN,4.0,4.0,
		1.0,0.0,NAN,2.0,NAN,1.0,3.0,NAN,5.0,5.0
	};

	List** L = getNeighbors_core(M, 0.5, 3, 2, S_FALSE, n, m);
	int** NIix = list2int(L, n);
	for (int j=0; j<n; j++)
		list_destroy(L[j]);
	free(L);

	// check neighbors of first row
	ASSERT_TRUE(NIix[0][0] == 3);
	for (int j=1; j<=3; j++)
	{
		ASSERT_TRUE(NIix[0][1] == j || 
			NIix[0][2] == j || NIix[0][3] == j);
	}
	for (int i=0; i<n; i++)
		free(NIix[i]);
	free(NIix);

	L = getNeighbors_core(M, -1.0, 3, 2, S_FALSE, n, m);
	NIix = list2int(L, n);
	for (int j=0; j<n; j++)
		list_destroy(L[j]);
	free(L);

	// check neighbors of fifth row
	ASSERT_TRUE(NIix[4][0] == 3);
	for (int j=2; j<=5; j++)
	{
		if (j == 4)
			continue; //not self-neighbor
		ASSERT_TRUE(NIix[4][1] == j || 
			NIix[4][2] == j || NIix[4][3] == j);
	}

	for (int i=0; i<n; i++)
		free(NIix[i]);
	free(NIix);
}
