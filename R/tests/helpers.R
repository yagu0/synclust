checkEquals = function(target, current, 
	tolerance=.Machine$double.eps^0.5, msg="", ...)
{
	#all.equal.numeric ?
	result = all.equal(target, current, tolerance=tolerance, ...)
	if (result != TRUE)
		stop(msg)
}

checkTrue = function(b, msg="")
{
	if (!b)
		stop(msg)
}
