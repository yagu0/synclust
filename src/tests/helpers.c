#include "tests/helpers.h"

//auxiliary to check vectors equality
int checkEqualV(double* v1, double* v2, int n)
{
	double epsilon = 1e-10; // arbitrary small value to make comparisons
	for (int i=0; i<n; i++)
	{
		if (fabs(v1[i] - v2[i]) >= epsilon)
			return S_FALSE;
	}
	return S_TRUE;
}

// auxiliary to count distinct values in an integer array
int countDistinctValues(int* v, int n)
{
	int maxVal = v[0];
	for (int i=1; i<n; i++)
	{
		if (v[i] > maxVal)
			maxVal = v[i];
	}
	int* kountArray = (int*)calloc(maxVal+1,sizeof(int));
	int res = 0;
	for (int i=0; i<n; i++)
	{
		if (kountArray[v[i]] == 0)
		{
			res++;
			kountArray[v[i]] = 1; // mark this value as "seen"
		}
	}
	free(kountArray);
	return res;
}

// helper to check clustering result
int labelIsProcessed(int label, int* processedLabels, int countProcessedLabels)
{
	for (int j=0; j<countProcessedLabels; j++)
	{
		if (processedLabels[j] == label)
			return S_TRUE;
	}
	return S_FALSE;
}

// check both clusters purity and size (ideally they all contain n/clustCount points)
int checkClustersProportions(int* clusters, int n, int clustCount, double tol)
{
	// initialize array to keep track of clusters labels
	int* processedLabels = (int*)malloc(clustCount*sizeof(int));
	for (int j=0; j<clustCount; j++)
		processedLabels[j] = -1;
	int countProcessedLabels = 0, clustSize = n/clustCount;

	int i=0;
	while (i<n)
	{
		// go to the next unprocessed label (if possible)
		while (i < n && labelIsProcessed(clusters[i], processedLabels, countProcessedLabels))
		{
			i++;
		}
		if (i >= n)
			break;

		int label = clusters[i];
		processedLabels[countProcessedLabels++] = label;

		// count elements in current cluster (represented by label)
		int countDataWithCurLabel = 0;
		for (int ii=0; ii<n; ii++)
		{
			if (clusters[ii] == label)
				countDataWithCurLabel++;
		}
		// if cluster is too small or too big, test fails
		if ( fabs(countDataWithCurLabel - clustSize) / n > tol )
		{
			free(processedLabels);
			return S_FALSE;
		}

		// now check counts per cluster (repartition);
		// the labels should not be spread between different (true) groups
		int maxCountLabelInClust = 0;
		for (int kounter=0; kounter<clustCount; kounter++)
		{
			int countLabelInClust = 0;
			for (int ii=kounter*clustSize; ii<(kounter+1)*clustSize; ii++)
			{
				if (clusters[ii] == label)
					countLabelInClust++;
			}
			if (countLabelInClust > maxCountLabelInClust)
				maxCountLabelInClust = countLabelInClust;
		}
		// we must have max(repartition) / clustSize >= 1 - tol
		if ((double)maxCountLabelInClust / clustSize < 1.0 - tol)
		{
			free(processedLabels);
			return S_FALSE;
		}
	}

	free(processedLabels);
	return S_TRUE;
}
