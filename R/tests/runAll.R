source("helpers.R")
source("t.clustering.R")
source("t.connexity.R")
source("t.utils.R")

dyn.load("../../src/synclust.so")

functions = c(lsf.str())
for (func in functions)
{
	#ou test length(grep("test.", function)) > 0
	if (nchar(func) > 5 && substr(func, 1, 5) == "test.")
	{
		print(paste("run",func))
		eval(call(func))
	}
}

#sample call for full package :
#t = findSyncVarRegions(method="convex",M=NULL,k=10,alpha=0.0,gmode=1,K=5,dtype="spath",cmeth="HC",pcoef=2.2,h=5e-4,eps=1e-5,maxit=3e3,showLL=TRUE,disp=TRUE)
