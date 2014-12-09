#include <R.h>
#include <Rdefines.h>
#include "sources/neighbors.h"
#include "sources/utils/algebra.h"
#include <cgds/List.h>

// Function to obtain neighborhoods.
// NOTE: alpha = weight parameter to compute distances; -1 means "adaptive"
// WARNING : M is given in columns
SEXP getNeighbors(
	SEXP M_, 
	SEXP k_, 
	SEXP alpha_, 
	SEXP gmode_, 
	SEXP simpleDists_
) {
	// get scalar arguments
	int k = INTEGER_VALUE(k_);
	double alpha = NUMERIC_VALUE(alpha_);
	int gmode = INTEGER_VALUE(gmode_);
	int simpleDists = LOGICAL_VALUE(simpleDists_);

	// extract infos from M and get associate pointer
	SEXP dim = getAttrib(M_, R_DimSymbol);
	int nrow = INTEGER(dim)[0];
	int ncol = INTEGER(dim)[1];
	// M is always given by columns: easier to process in rows
	double* pM = transpose(REAL(M_), nrow, ncol);

	// Main call to core algorithm which fills neighborhoods lists
	List** neighborhoods = getNeighbors_core(pM, alpha, k, gmode, simpleDists, nrow, ncol);

	// transfer neighborhoods lists into R vectors
	SEXP NIix, NIds;
	PROTECT(NIix = allocVector(VECSXP, nrow)); //indices
	PROTECT(NIds = allocVector(VECSXP, nrow)); //distances
	for (int i=0; i<nrow; i++)
	{
		SEXP neighbsIX, neighbsDS;
		PROTECT(neighbsIX = NEW_INTEGER(list_size(neighborhoods[i])));
		PROTECT(neighbsDS = NEW_NUMERIC(list_size(neighborhoods[i])));
		int* pNeighbsIX = INTEGER_POINTER(neighbsIX);
		double* pNeighbsDS = NUMERIC_POINTER(neighbsDS);
		ListIterator* neighbsI = list_get_iterator(neighborhoods[i]);
		int j = 0;
		while (listI_has_data(neighbsI))
		{
			IndDist indDist; listI_get(neighbsI, indDist);
			// WARNING: R arrays start at index 1
			pNeighbsIX[j] = indDist.index + 1;
			pNeighbsDS[j] = indDist.distance;
			j++;
			listI_move_next(neighbsI);
		}
		SET_VECTOR_ELT(NIix, i, neighbsIX);
		SET_VECTOR_ELT(NIds, i, neighbsDS);
		UNPROTECT(2);
		listI_destroy(neighbsI);
		list_destroy(neighborhoods[i]);
	}
	free(neighborhoods);

	// create R list labels to access with NI$ix and NI$ds
	SEXP listNames;
	char* lnames[2] = {"ix", "ds"}; //lists labels
	PROTECT(listNames = allocVector(STRSXP,2));
	for (int i=0; i<2; i++)
		SET_STRING_ELT(listNames,i,mkChar(lnames[i]));

	// allocate and fill neighborhoods list to return
	SEXP NI;
	PROTECT(NI = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(NI, 0, NIix);
	SET_VECTOR_ELT(NI, 1, NIds);
	setAttrib(NI, R_NamesSymbol, listNames);

	UNPROTECT(4);
	return NI;
}
