#test matrix [de]standardization
test.de_standardize = function()
{
	n = 100
	m = 10
	M = matrix(rnorm(n*m,mean=2.0,sd=2.0),nrow=n,ncol=m)
	
	std = synclust:::standardize(M)
	#result is centered:
	checkEquals(rep(0.0, m), colMeans(std$M))
	#result is standardized:
	checkEquals(rep(1.0, m), sqrt( unlist( apply(std$M, 2, var) ) ))
	
	#rebuilt M == M:
	M_rec = synclust:::destandardize(std)
	checkEquals(M, M_rec)
}

#test neighborhoods remapping into one smaller component
test.remap_neighbors = function()
{
	#connex comp : 1-2-5-8-10
	#to remap into 1-2-3-4-5
	NI = list(
		"ix" = list(
			c(2,8,10), #V(1)
			c(1,5,8,10), #V(2)
			c(4,7,11), #V(3)
			c(3,6,9,12), #V(4)
			c(2,10), #V(5)
			c(4,7), #V(6)
			c(3,6,9,12), #V(7)
			c(1,2,10), #V(8)
			c(4,7,11), #V(9)
			c(1,2,5,8), #V(10)
			c(3,9), #V(11)
			c(4,7)), #V(12)
		"ds" = list(
			c(1.0,2.0,3.0), #1
			c(1.0,2.0,3.0,4.0), #2
			c(1.0,1.0,1.0), #3
			c(1.0,1.0,1.0,1.0), #4
			c(2.0,2.0), #5
			c(1.0,1.0), #6
			c(1.0,1.0,1.0,1.0), #7
			c(2.0,3.0,1.0), #8
			c(1.0,1.0,1.0), #9
			c(3.0,4.0,2.0,1.0), #10
			c(1.0,1.0), #11
			c(1.0,1.0))) #12

	indices = c(1,2,5,8,10)
	locNI = synclust:::remapNeighbors(NI, indices)
	checkEquals(2, length(locNI))
	checkEquals(length(indices), length(locNI$ix))
	checkEquals(length(indices), length(locNI$ds))
	
	#formerly index 1 (now 1)
	checkEquals(c(2,4,5), locNI$ix[[1]])
	checkEquals(NI$ds[[1]], locNI$ds[[1]],)
	#formerly index 2 (now 2)
	checkEquals(c(1,3,4,5), locNI$ix[[2]])
	checkEquals(NI$ds[[2]], locNI$ds[[2]])
	#formerly index 5 (now 3)
	checkEquals(c(2,5), locNI$ix[[3]])
	checkEquals(NI$ds[[5]], locNI$ds[[3]])
	#formerly index 8 (now 4)
	checkEquals(c(1,2,5), locNI$ix[[4]])
	checkEquals(NI$ds[[8]], locNI$ds[[4]])
	#formerly index 10 (now 5)
	checkEquals(c(1,2,3,4), locNI$ix[[5]])
	checkEquals(NI$ds[[10]], locNI$ds[[5]])
}
