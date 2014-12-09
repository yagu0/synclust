#bipartite graph
test.connexity1 = function()
{
	n = 10
	NIix = as.list(rep(NA,n))
	#connect 0 with 1, 2 with 3 ...
	for (i in 2*(0:(n/2-1)) + 1)
	{
 		NIix[[i]] = i+1
 		NIix[[i+1]] = i
	}
	cc = synclust:::getConnectedComponents(NIix)
	#cc should contain exactly n/2 integers
	checkEquals(n/2, length(unique(cc)))
}

#cyclic graph
test.connexity2 = function()
{
	n = 10
	NIix = as.list(rep(NA,n))
	#connect 0 with 1, 1 with 2 ...
	for (i in 1:n)
 		NIix[[i]] = c(ifelse(i==1,n,i-1), i%%n+1)
	cc = synclust:::getConnectedComponents(NIix)
	#cc should contain only one integer (1)
	checkEquals(1, length(unique(cc)))
}

#custom graph with 3 connex components
test.connexity3 = function()
{
	n = 10
	NIix = as.list(rep(0,n))
	NIix[[1]] = c(3,5)
	NIix[[2]] = c(3,5)
	NIix[[3]] = c(1,2)
	NIix[[4]] = c(6,9,10)
	NIix[[5]] = c(1,2)
	NIix[[6]] = c(4)
	NIix[[7]] = c(8)
	NIix[[8]] = c(7)
	NIix[[9]] = c(4)
	NIix[[10]] = c(4,9)
	cc = synclust:::getConnectedComponents(NIix)
	#cc should contain only three integers
	checkEquals(3, length(unique(cc)))
}

#custom graph, 1 connex component
test.connexity4 = function()
{
	n = 10
	NIix = as.list(rep(0,n))
	NIix[[1]] = c(3,4,8)
	NIix[[2]] = c(3,5,7)
	NIix[[3]] = c(1,2)
	NIix[[4]] = c(1,6,9,10)
	NIix[[5]] = c(2)
	NIix[[6]] = c(4,8)
	NIix[[7]] = c(2)
	NIix[[8]] = c(1,6,10)
	NIix[[9]] = c(4)
	NIix[[10]] = c(4,8)
	cc = synclust:::getConnectedComponents(NIix)
	#cc should contain only one integer (1)
	checkEquals(1, length(unique(cc)))
}
