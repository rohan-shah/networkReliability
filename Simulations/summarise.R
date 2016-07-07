source("./generateScenarios.R")
library(Rmpfr)
library(networkReliability)
allResults <- list()
for(i in 1:nrow(scenarios))
{
	path <- file.path("results", scenarios[i, "file"])
	if(file.exists(path))
	{
		load(path)
		allResults[[i]] <- results
	}
	else allResults[[i]] <- list()
}
secondsPerRun <- lapply(allResults, function(x) unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs"))))
averageSecondsPerRun <- unlist(lapply(secondsPerRun, function(x) if(is.null(x)) NA else mean(x)))
averageEstimatesFunc <- function(x)
{
	if(length(x) == 0)
	{
		return(NA)
	}
	else if(class(x[[1]]) == "approximateZeroVarianceResult") 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@firstMoment)))))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else 
	{
		stop("Internal error")
	}
}
averageEstimates <- do.call(c, lapply(allResults, averageEstimatesFunc))

varianceFunc <- function(x)
{
	if(length(x) == 0)
	{
		return(NA)
	}
	else if(class(x[[1]]) == "approximateZeroVarianceResult") 
	{
		return(var(as.numeric(do.call(c, lapply(x, function(y) y@firstMoment)))))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORResult")
	{
		return(var(as.numeric(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else 
	{
		stop("Internal error")
	}
}
variances <- do.call(c, lapply(allResults, varianceFunc))
workNormalizedVariance <- as.numeric(variances * averageSecondsPerRun)

data(dodecahedronConnectedEnumeration)
trueValues <- apply(scenarios, 1, function(x)
{
	if(x["graph"] == "dodecahedron")
	{
		return(exhaustiveProbability(searchObj = dodecahedronConnectedEnumeration, probability = as.numeric(x["probability"])))
	}
	return(NA)
})
empiricalBias <- averageEstimates - trueValues
relativeErrors <- sqrt(variances) / averageEstimates
wnrv <- as.numeric(variances * averageSecondsPerRun / (averageEstimates^2))

save(allResults,averageEstimates, averageSecondsPerRun, empiricalBias, secondsPerRun, variances, workNormalizedVariance, wnrv, relativeErrors, trueValues, file = "summarised.RData")
