#build similarity matrix W (NOTE : sparse matrix ==> optimizations later)
getSimilarityMatrix = function(NI)
{
	# using a local sigma would be nice, but break W symmetry,
	# which cannot easily be repaired then (??!)
	# ==> we use a global sigma, with a very simple heuristic
	
	n = length(NI$ix)
	distances = c()
	for (i in 1:n) distances = c(distances,NI$ds[[i]])
	distances = unique(distances)
	sigma2 = median(distances)^2 #for example...

	W = matrix(0.0,nrow=n,ncol=n)
	for (i in 1:n)
		W[ i, NI$ix[[i]] ] = exp( - NI$ds[[i]]^2 / sigma2 )

	return (W)
}

#epsilon constant, used as a zero threshold
EPS = 100 * .Machine$double.eps

#Moore-Penrose pseudo inverse
mppsinv = function(M)
{
	s = svd(M)
	sdiag = s$d ; sdiag[sdiag < EPS] = Inf
	p = min(nrow(M),ncol(M))
	sdiag = diag(1.0 / sdiag, p)
	return ((s$v) %*% sdiag %*% t(s$u))
}

#get distance matrix from data and similarity : Commute Time
getECTDistances = function(NI)
{
	n = length(NI$ix) ; seqVect = 1:n
	if (n <= 1) return (0.0) #nothing to do...
	
	#get laplacian (...inverse) :
	W = getSimilarityMatrix(NI)
	invLap = mppsinv(diag(rowSums(W)) - W)

	#...and distances
	ectd = matrix(0.0, nrow=n, ncol=n)
	for (ij in 1:n)
	{
		ectd[ij,] = ectd[ij,] + invLap[ij,ij]
		ectd[,ij] = ectd[,ij] + invLap[ij,ij]
	}
	ectd = ectd - 2*invLap
	return (ectd)
}

# Call Dijsktra algorithm on every vertex to build distances matrix
getShortestPathDistances = function(NI)
{
	n = length(NI$ix)
	distancesIn = matrix(NA,nrow=n,ncol=n)
	for (i in 1:n)
		distancesIn[i,NI$ix[[i]]] = NI$ds[[i]]

	distancesOut = matrix(nrow=n, ncol=n)
	for (i in 1:n)
		distancesOut[i,] = .Call("dijkstra", distancesIn, i)
	return (distancesOut)
}

## MAIN CALL to get distances matrix
getDistances = function(dtype, NI)
{
	distances = matrix()
	if (dtype=="spath")
		distances = getShortestPathDistances(NI)
	else if (dtype=="ectd")
		distances = getECTDistances(NI)
	
	diag(distances) = 0.0 #distances to self are zero
	return (distances)
}
