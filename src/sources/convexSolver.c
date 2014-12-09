#include "sources/convexSolver.h"
#include <stdio.h> //to trace LL evolution
#include <stdlib.h>
#include <math.h>
#include "sources/utils/algebra.h"

// auxiliary to compute log-likelihood + penalty
double computeLogLikelihood(
	double** f, double* theta, double** Z, double*** phi, 
	int* lengthNIix, int** NIix, double alpha, int nrow, int ncol) 
{
	double LL = 0.0;

	// for each row in data matrix: (one row = observations from 2001 to 2009)
	for (int i=0; i<nrow; i++)
	{
		// theta[i] == -INFINITY if no birds were seen at this site
		if (theta[i] != -INFINITY)
		{
			// for each year from 2001 to 2009:
			for (int j=0; j<ncol; j++)
				LL += (exp(theta[i] + f[i][j]) - Z[i][j] * (theta[i] + f[i][j]));
		}
		// add penalty term
		double penalty = 0.0;
		double Ds = 1.0/lengthNIix[i];
		// lengthNIix[i] == size of the neighborhood of site i
		for (int j=0; j<lengthNIix[i]; j++)
		{
			// compute <phi[s,u] , f[s,] - f[u,]> with u == NIix[i][j] (j-th neighbor of i)
			double dotProduct = 0.0;
			for (int jj=0; jj<ncol; jj++)
				dotProduct += phi[i][NIix[i][j]][jj] * (f[i][jj] - f[NIix[i][j]][jj]);
			// normalization by sum of inverses of neighborhoods sizes
			double Dsu = Ds + 1.0/lengthNIix[NIix[i][j]];
			penalty += dotProduct / Dsu;
		}
		LL += alpha * penalty;
	}

	return LL;
}

// compute estimated ("repaired", "smoothed"...) variations from rows of M
// NOTE: geographic coordinates dropped here, since they are unused
Parameters getVarsWithConvexOptim_core(
			double* pM, int* lengthNIix, int** NIix, int nrow, int ncol, 
			double alpha, double h, double epsilon, int maxiter, bool symmNeighbs, bool trace)
{
	double EPS = 1e-10; // HACK: some small numerical constant to avoid oddities

	// theta_s = log(average z_st)
	double* theta = (double*)calloc(nrow,sizeof(double));
	// NOTE:: Z is 'double' because it is [can be] an average value (of integers)
	double** Z = (double**)malloc(nrow*sizeof(double*));
	for (int i=0; i<nrow; i++)
	{
		Z[i] = (double*)malloc(ncol*sizeof(double));
		for (int j=0; j<ncol; j++)
		{
			Z[i][j] = pM[i*ncol+j];
			theta[i] += Z[i][j];
		}
		// since pM values are assumed to be integers (and ncol not too high ?!),
		// the following test may be simplified into (theta[i]==0.0)
		if (fabs(theta[i]) < EPS)
			theta[i] = -INFINITY;
		else
			theta[i] = log(theta[i]/ncol);
	}
	// initialize f to observed variations
	double** F = (double**)malloc(nrow*sizeof(double*));
	for (int i=0; i<nrow; i++)
	{
		F[i] = (double*)calloc(ncol,sizeof(double));
		if (theta[i] != -INFINITY)
		{
			for (int j=0; j<ncol; j++)
			{
				if (Z[i][j] > 0.0)
					F[i][j] = log(Z[i][j]) - theta[i];
			}
		}
	}
	// phi_s,u = 1/sqrt(ncol) (1 ... 1) [TODO:: costly in memory !]
	double invSqrtNcol = 1.0/sqrt(ncol);
	double*** phi = (double***)malloc(nrow*sizeof(double**));
	for (int i=0; i<nrow; i++)
	{
		phi[i] = (double**)malloc(nrow*sizeof(double*));
		for (int ii=0; ii<nrow; ii++)
		{
			phi[i][ii] = (double*)malloc(ncol*sizeof(double));
			for (int j=0; j<ncol; j++)
				phi[i][ii][j] = invSqrtNcol;
		}
	}

	// initialize log-likelihood
	double LL = computeLogLikelihood(
		F, theta, Z, phi, lengthNIix, NIix, alpha, nrow, ncol);
	double oldLL = -INFINITY;

	/*******************
	 * START ITERATIONS
	 *******************/

	int kounter = 0; // limit iterations count, in case of
	while (fabs(LL - oldLL) >= epsilon && kounter++ < maxiter)
	{
		// gradient descent for theta
		for (int i=0; i<nrow; i++) {
			if (theta[i] == -INFINITY)
				continue; // skip these sites: cannot get information
			double sumExpFst = 0.0;
			for (int j=0; j<ncol; j++)
				sumExpFst += exp(F[i][j]);
			double sumZst = 0.0;
			for (int j=0; j<ncol; j++)
				sumZst += Z[i][j];
			double gradI = exp(theta[i]) * sumExpFst - sumZst;
			theta[i] -= h * gradI;
		}

		// gradient descent for f
		double sumDdivPhi;
		for (int i=0; i<nrow; i++)
		{
			double invDs = 1.0/lengthNIix[i];
			for (int j=0; j<ncol; j++)
			{
				double gradIJ = - Z[i][j];
				if (theta[i] != -INFINITY)
				{
					// no theta[i] contribution if nothing observed
					gradIJ += exp(theta[i] + F[i][j]);
				}
				// + sum on neighbors u: s-->u, - sum on neighbors u: u-->s
				sumDdivPhi = 0.0;
				for (int jj=0; jj<lengthNIix[i]; jj++)
				{
					double Dsu = invDs + 1.0/lengthNIix[NIix[i][jj]];
					sumDdivPhi += phi[i][NIix[i][jj]][j] / Dsu;
					if (symmNeighbs)
					{
						//shortcut: if symmetric neighborhoods, it's easy to sum on u-->s
						sumDdivPhi -= phi[NIix[i][jj]][i][j] / Dsu;
					}
				}
				gradIJ += alpha * sumDdivPhi;
				if (!symmNeighbs)
				{
					// need to remove sum on neighbors u: u-->s, uneasy way.
					//TODO: computation is much too costly here; need preprocessing
					sumDdivPhi = 0.0;
					for (int ii=0; ii<nrow; ii++)
					{
						//~ if (ii == i) continue;
						for (int jj=0; jj<lengthNIix[ii]; jj++)
						{
							if (NIix[ii][jj] == i)
							{
								sumDdivPhi += phi[ii][i][j] / (invDs + 1.0/lengthNIix[ii]);
								break; //i can only appear once among neighbors of ii
							}
						}
					}
					gradIJ -= alpha * sumDdivPhi;
				}
				F[i][j] -= h * gradIJ;
			}
		}

		// gradient ascent for phi
		for (int i=0; i<nrow; i++)
		{
			double Ds = 1.0/lengthNIix[i];
			for (int ii=0; ii<nrow; ii++)
			{
				double Dsu = Ds + 1.0/lengthNIix[ii];
				for (int j=0; j<ncol; j++)
				{
					double gradI_II_J = alpha * (F[i][j] - F[ii][j]) / Dsu;
					phi[i][ii][j] += h * gradI_II_J;
				}
				// also renormalize to have ||phi_su|| == 1.0
				double normPhi = norm2(phi[i][ii], ncol);
				//~ if (normPhi > 1.0) {
				if (normPhi > EPS)
				{
					for (int j=0; j<ncol; j++)
						phi[i][ii][j] /= normPhi;
				}
			}
		}

		oldLL = LL;
		LL = computeLogLikelihood(
			F, theta, Z, phi, lengthNIix, NIix, alpha, nrow, ncol);
		if (trace)
			printf("%i / LLs: %.0f %.0f\n",kounter,oldLL,LL); // optional trace of LL evolution
	} /*** END ITERATIONS ***/

	// free all local parameters arrays but (theta, F) (used as return value)
	for (int i=0; i<nrow; i++)
	{
		free(Z[i]);
		for (int ii=0; ii<nrow; ii++)
			free(phi[i][ii]);
		free(phi[i]);
	}
	free(Z);
	free(phi);

	Parameters params;
	params.f = F;
	params.theta = theta;
	return params;
}
