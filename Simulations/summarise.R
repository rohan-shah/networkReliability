source("./generateScenarios.R")
library(Rmpfr)
allResults <- list()
for(i in 1:nrow(scenarios))
{
	path <- file.path("results", scenarios[i, "file"])
	if(file.exists(path))
	{
		load(path)
		allResults[[i]] <- results
	}
}
secondsPerRun <- lapply(allResults, function(x) unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs"))))
averageSecondsPerRun <- unlist(lapply(secondsPerRun, function(x) if(is.null(x)) NA else mean(x)))
secondsSingleSample <- unlist(lapply(allResults, function(x)
	{
		if(is.null(x) || class(x[[1]]) == "generalisedSplittingResult") 
		{
			return(NA)
		}
		else 
		{
			return(mean(unlist(lapply(x, function(y) difftime(y@end, y@start, units = "secs")/y@n))))
		}
	}))

averageEstimatesFunc <- function(x)
{
	if(is.null(x))
	{
		return(NA)
	}
	if(class(x[[1]]) == "crudeMCResult") 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@data)))))
	}
	else if(class(x[[1]]) %in% c("generalisedSplittingFixedFactorsResult", "generalisedSplittingFixedEffortResult"))
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@firstMomentSingleSample)))))
	}
}
averageEstimates <- do.call(c, lapply(allResults, averageEstimatesFunc))

varianceSingleSampleFunc <- function(x)
{
	if(is.null(x))
	{
		return(NA)
	}
	if(class(x[[1]]) == "crudeMCResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@data*(1 - y@data))))))
	}
	else if(class(x[[1]]) == "generalisedSplittingFixedFactorsResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@estimatedVariance)))) / x[[1]]@n)
	}
	else if(class(x[[1]]) == "generalisedSplittingFixedEffortResult")
	{
		return(NA)
	}
	else 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@varianceSingleSample)))))
	}
}
varianceSingleSample <- do.call(c, lapply(allResults, varianceSingleSampleFunc))
workNormalizedSingleSampleVariance <- as.numeric(varianceSingleSample * secondsSingleSample)

varianceFunc <- function(x)
{
	if(is.null(x))
	{
		return(NA)
	}
	if(class(x[[1]]) == "crudeMCResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@data*(1 - y@data)/y@n)))))
	}
	else if(class(x[[1]]) == "generalisedSplittingFixedFactorsResult")
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@estimatedVariance)))))
	}
	else if(class(x[[1]]) ==  "generalisedSplittingFixedEffortResult")
	{
		return(as.numeric(var(do.call(c, lapply(x, function(y) y@estimate)))))
	}
	else 
	{
		return(as.numeric(mean(do.call(c, lapply(x, function(y) y@varianceSingleSample/y@n)))))
	}
}
variances <- do.call(c, lapply(allResults, varianceFunc))
workNormalizedVariance <- as.numeric(variances * averageSecondsPerRun)

relativeErrors <- sqrt(variances) / averageEstimates
wnrv <- as.numeric(variances * averageSecondsPerRun / (averageEstimates^2))

save(allResults, secondsSingleSample, varianceSingleSample, workNormalizedSingleSampleVariance, averageEstimates, averageSecondsPerRun, secondsPerRun, variances, workNormalizedVariance, wnrv, relativeErrors, file = "summarised.RData")
