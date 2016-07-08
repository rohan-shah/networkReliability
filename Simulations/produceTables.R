source("./generateScenarios.R")
library(xtable)
load("./summarised.RData")
produceTable <- function(indices, showEmpiricalMse)
{
	mse <- as.numeric(mse)
	results <- cbind(scenarios, wnrv, relativeErrors, mse)
	if(showEmpiricalMse)
	{
		retVal <- results[indices,c("method", "sampleSize", "wnrv", "relativeErrors", "mse")]
		colnames(retVal) <- c("Method", "$n$", "WNRV", "RE", "MSE")
	}
	else
	{
		retVal <- results[indices,c("method", "sampleSize", "wnrv", "relativeErrors")]
		colnames(retVal) <- c("Method", "$n$", "WNRV", "RE")
	}
	retVal$Method <- c("approxZeroVarianceWOR" = "WOR", "approxZeroVariance" = "IS")[retVal$Method]
	return(retVal)
}
indices <- which(scenarios$graph == "dodecahedron" & scenarios$probability == 0.99)
print(xtable(produceTable(indices, TRUE), display = c("s", "s", "d", "e", "e", "e"), label = "table:networkReliability_dodec", caption = "Simulation results for the dodecahedron graph with $p = 0.99$."), include.rownames=FALSE, math.style.exponents = TRUE, sanitize.text.function = identity)

indices <- which(scenarios$graph == "grid9Augmented" & scenarios$probability == 0.9999)
print(xtable(produceTable(indices, FALSE), display = c("s", "s", "d", "e", "e"), label = "table:networkReliability_grid9Augmented", caption = "Simulation results for the augmented $9\\times 9$ grid graph with $p = 0.9999$."), include.rownames=FALSE, math.style.exponents = TRUE, sanitize.text.function = identity)
