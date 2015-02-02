#' Direct clustering from a neighborhoods graph, or get regions from (Poisson) 
#' distribution parameters optimization, using convex relaxation.
#'
#' @param method Global method: "direct" or "convex"
#' @param M Matrix of observations in rows, the two last columns 
#'        corresponding to geographic coordinates; 
#'        set to NULL to use our initial dataset (625 rows / 9 years)
#' @param k Number of neighbors
#' @param alpha Weight parameter for intra-neighborhoods distance computations; 
#'        0 = take only geographic coordinates into account; 
#'        1 = take only observations over the years into account; 
#'        in-between : several levels of compromise; 
#'        -1 or any negative value : use a heuristic to choose alpha.
#' @param gmode Neighborhood type. 0 = reduced [mutual] kNN; 1 = augmented kNN (symmetric); 
#'        2 = normal kNN; 3 = one NN in each quadrant; (NON-symmetric). 
#'        NOTE: gmode==3 automatically sets k==4 (at most!)
#' @param K Number of clusters
#' @param dtype Distance type, in {"simple","spath","ectd"}. 
#'        NOTE: better avoid "simple" if gmode>=2
#' @param cmeth Clustering method, in {"KM","HC","spec"} for k-means (distances based) 
#'        or hierarchical clustering, or spectral clustering (only if gmode>=2)
#' @param pcoef Penalty value for convex optimization [default: 1.0]
#' @param h Step in the min LL algorithm [default: 1e-3]
#' @param eps Threshold to stop min.LL iterations [default: 1e-3]
#' @param maxit Maximum number of iterations in the min LL algo [default: 1e3]
#' @param showLL Print trace of log-likelihood evolution [default: true]
#' @param disp True [default] for interactive display (otherwise nothing gets plotted)
#' @return list with the following entries. M: data matrix in input; NI: computed neighborhoods; 
#'         dists: computed distances matrix; clusts: partition into K clusters, as an integer vector; 
#'         cxpar: parameters obtained after convex optimization (if applicable)
#' @export
#' @examples
#' cvr = findSyncVarRegions("convex",M=NULL,k=10,alpha=0.1,gmode=1,K=5,dtype="spath",cmeth="HC")
#' drawMapWithSitez(cvr$M, cvr$clusters)
#' drawNeighboroodGraph(cvr$M, cvr$NI)
#'
findSyncVarRegions = function(method, M, k, alpha, gmode, K, dtype, cmeth,
	pcoef=1.0, h=1e-3, eps=1e-3, maxit=1e3, showLL=TRUE, disp=TRUE)
{
	#get matrix M if not directly provided
	if (is.null(M))
	{
		data("example", package="synclust")
		M = synclust_sample
	}
	if (is.character(M))
		M = as.matrix(read.table(M))

	n = nrow(M)
	m = ncol(M)

	#pretreatment for neighborhoods search: standardize M columns
	#TODO: maybe apply only on coordinates columns ?
	std = standardize(M)

	#get neighborhoods [FALSE because NOT simpleDists; see C code]
	NI = .Call("getNeighbors", std$M, k, alpha, gmode, FALSE)

	#optional intermediate display : map + graph (monocolor)
	if (disp)
		promptForMapDisplay("interm", M[,(m-1):m], NIix=NI$ix)

	clusters = rep(1,n)
	distances = matrix(NA,nrow=n,ncol=n)
	cxpar = list()

	## DIRECT CLUSTERING ##
	if (method=="direct")
	{
		if (gmode >= 2)
			stop("'gmode' must be 0 or 1 for direct clustering")
		if (dtype=="simple")
			stop("'dtype' cannot be set to \"simple\" for direct (graph!) clustering")

		#find connected components in the graph defined by NI
		cc = reordering(.Call("getConnectedComponents", NI$ix))
		nbC = max(cc)
		if (nbC > 10)
			stop(paste("ABORT: too many connex components (found ",nbC,")",sep=''))
		if (nbC > 1)
			print(paste("*** WARNING:",nbC,"connex components ***"))
		clusters = cc

		#for each connected component...
		for (i in 1:nbC)
		{
			indices = (1:n)[cc == i]
			nc = length(indices)
			if (nc <= 1)
				next #nothing to do with such a small component
			
			if (nbC > 1)
			{
				doClust = readline(paste(">>> cluster current component of cardinal",nc,"/",n,"? (y/n)\n"))
				if (doClust == "y")
					K = readline(">>> into how many groups ? (int >= 2)\n")
				else
					next
			}
			
			#0] remap NI in current connex component
			locNI = remapNeighbors(NI, indices)
			
			#1] determine similarity and distance matrices (e.g. using a random walk)
			locDists = getDistances(dtype, locNI)
			distances[indices,indices] = locDists
			
			#2] cluster data inside connex component according to distance matrix
			locClusters = getClusters(locDists, cmeth, K)
			maxInd = max(clusters)
			clusters[indices] = locClusters + maxInd #avoid indices overlaps
		}
	}

	## CONVEX RELAXATION ##
	else if (method=="convex")
	{		
		#preliminary: remove NA's by averaging over each serie's values
		M = replaceNAs(M)
		
		#use NI$ix and matrix M to estimate initial parameters,
		#and then iterate until convergence to get f + theta
		#theta == mean observations count at each site s
		#f == estimated variations at each site ("time-series" of T points)
		cxpar = .Call("getVarsWithConvexOptim", 
			M[,1:(m-2)], NI$ix, pcoef, h, eps, maxit, (gmode <= 1), showLL)
		f = cxpar$f #the only one we use (others can be checked by user later)
		
		#cluster "time-series" f, using simple kmeans/HC, spect.clust, 
		#or [in a graph] KM or HC, after redefining a NI (using f only)
		
		if (dtype=="simple")
		{
			#use R core functions
			if (cmeth=="KM")
				clusters = kmeans(f, K, iter.max=100, nstart=10)$cluster
			else if (cmeth=="HC")
			{
				hct = hclust(dist(f), method="ward")
				clusters = cutree(hct, K)
			}
			else if (cmeth=="spec")
			{
				require(kernlab)
				clusters = as.integer(specc(f, K, kpar="automatic"))
			}
		}
		
		else
		{
			# recompute NI from repaired/smoothed data [simpleDists=TRUE, no graph dists]
			#NOTE: gmode==1, augmented kNN (arbitrary, but should be symmetric)
			NI = .Call("getNeighbors", f, k, alpha, 1, TRUE)
			
			#find connected components in the graph defined by NI
			cc = reordering(.Call("getConnectedComponents", NI$ix))
			
			nbC = max(cc)
			if (nbC > 10) 
				stop(paste("ABORT: too many connex components (found ",nbC,")",sep=''))
			if (nbC > 1) 
				print(paste("*** WARNING:",nbC,"connex components ***"))
			clusters = cc
			
			#for each connected component...
			for (i in 1:nbC)
			{
				indices = (1:n)[cc == i]
				nc = length(indices)
				if (nc <= 1)
					next #nothing to do with such a small component
				
				if (nbC > 1)
				{
					doClust = readline(paste(">>> cluster current component of cardinal",nc,"/",n,"? (y/n)\n"))
					if (doClust == "y")
						K = readline(">>> into how many groups ? (int >= 2)\n")
					else
						next
				}
				
				#0] remap NI in current connex component
				locNI = remapNeighbors(NI, indices)
				
				#1] determine similarity and distance matrices (e.g. using a random walk)
				locDists = getDistances(dtype, locNI)
				distances[indices,indices] = locDists
				
				#2] cluster data inside connex component according to distance matrix
				locClusters = getClusters(locDists, cmeth, K)
				maxInd = max(clusters)
				clusters[indices] = locClusters + maxInd #avoid indices overlaps
			}
		}
	}

	clusters = reordering(clusters)
	#optional final display : map with clusters colors
	if (disp)
		promptForMapDisplay("final", M[,(m-1):m], clusters=clusters)

	#give back matrix M as given to the function
	M = destandardize(std)

	return (list("M"=M, "NI"=NI, "dists"=distances, "clusts"=clusters, "cxpar"=cxpar))
}
