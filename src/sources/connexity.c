#include "sources/connexity.h"
#include <cgds/Stack.h>
#include <stdlib.h>
#include "sources/utils/boolean.h"

int* getConnectedComponents_core(int** NIix, int* lengthNIix, int n)
{
	int* cc = (int*)calloc(n,sizeof(int));
	Stack* toBeExplored = stack_new(int);
	int curInd = 0, nextInd;
	bool* alreadyExpanded = (bool*)calloc(n,sizeof(bool));

	// while the entire graph hasn't been visited
	while (S_TRUE)
	{
		int label = curInd+1;
		cc[curInd] = label;
		stack_push(toBeExplored, curInd);

		// while new elements are discovered in current component,
		// mark them as expanded and stack their neighbors
		while (!stack_empty(toBeExplored))
		{
			stack_top(toBeExplored, nextInd);
			stack_pop(toBeExplored);
			cc[nextInd] = label;

			for (int j=0; j<lengthNIix[nextInd]; j++)
			{
				if (!alreadyExpanded[NIix[nextInd][j]])
					stack_push(toBeExplored, NIix[nextInd][j]);
			}
			alreadyExpanded[nextInd] = S_TRUE;
		}

		// curInd is set to the next unexplored index (if possible)
		for (int i=0; i<n; i++)
		{
			if (cc[i] == 0)
			{
				curInd = i;
				break;
			}
		}
		if (cc[curInd] != 0)
			break;
	}

	free(alreadyExpanded);
	stack_destroy(toBeExplored);

	return cc;
}
