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
	else if(class(x[[1]]) == "residualResamplingResult") 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@firstMoment)))))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORWithVarianceResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else 
	{
		stop("Internal error")
	}
}
averageEstimates <- do.call(c, lapply(allResults, averageEstimatesFunc))

singleEstimatesFunc <- function(x)
{
	if(length(x) == 0)
	{
		return(NA)
	}
	else if(class(x[[1]]) == "residualResamplingResult") 
	{
		return(as.numeric(x[[1]]@firstMoment))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceResult") 
	{
		return(as.numeric(x[[1]]@firstMoment))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORResult")
	{
		return(as.numeric(x[[1]]@estimate))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORWithVarianceResult")
	{
		return(as.numeric(x[[1]]@estimate))
	}
	else 
	{
		stop("Internal error")
	}
}
singleEstimates <- do.call(c, lapply(allResults, singleEstimatesFunc))

empiricalVarianceFunc <- function(x)
{
	if(length(x) == 0)
	{
		return(NA)
	}
	else if(class(x[[1]]) == "residualResamplingResult") 
	{
		return(var(as.numeric(do.call(c, lapply(x, function(y) y@firstMoment)))))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceResult") 
	{
		return(var(as.numeric(do.call(c, lapply(x, function(y) y@firstMoment)))))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORResult")
	{
		estimates <- do.call(c, lapply(x, function(y) y@estimate))
		return(as.numeric(sum(estimates*estimates)/length(estimates) - (sum(estimates)/length(estimates))^2))
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORWithVarianceResult")
	{
		estimates <- do.call(c, lapply(x, function(y) y@estimate))
		return(as.numeric(sum(estimates*estimates)/length(estimates) - (sum(estimates)/length(estimates))^2))
		#return(var(as.numeric(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else 
	{
		stop("Internal error")
	}
}
empiricalVariances <- do.call(c, lapply(allResults, empiricalVarianceFunc))

varianceEstimateFunc <- function(x)
{
	if(length(x) == 0)
	{
		return(NA)
	}
	else if(class(x[[1]]) == "residualResamplingResult") 
	{
		return(NA)
	}
	else if(class(x[[1]]) == "approximateZeroVarianceResult") 
	{
		return(NA)
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORResult")
	{
		return(NA)
	}
	else if(class(x[[1]]) == "approximateZeroVarianceWORWithVarianceResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@varianceEstimate)))))
	}
	else 
	{
		stop("Internal error")
	}
}
varianceEstimates <- do.call(c, lapply(allResults, varianceEstimateFunc))
workNormalizedVariance <- as.numeric(empiricalVariances * averageSecondsPerRun)
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
relativeErrors <- sqrt(empiricalVariances) / averageEstimates
wnrv <- as.numeric(empiricalVariances * averageSecondsPerRun / (averageEstimates^2))

mse <- empiricalVariances + empiricalBias^2
save(averageEstimates, singleEstimates, varianceEstimates, mse,  averageSecondsPerRun, empiricalBias, secondsPerRun, empiricalVariances, workNormalizedVariance, wnrv, relativeErrors, trueValues, file = "summarised.RData")
