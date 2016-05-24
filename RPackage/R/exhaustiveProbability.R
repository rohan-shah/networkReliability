exhaustiveProbability <- function(searchObj, probability)
{
	if(!isS4(searchObj) || !is(searchObj, "exhaustiveSearchResult"))
	{
		stop("Input searchObj must be an object of class exhaustiveSearchResult")
	}
	.Call("exhaustiveProbability", searchObj, probability, PACKAGE="networkReliability")
}
