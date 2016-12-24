source("./generateScenarios.R")
library(xtable)
load("./summarised.RData")
produceTable <- function(indices, showEmpiricalMse)
{
	mse <- as.numeric(mse)
	results <- cbind(scenarios, wnrv, relativeErrors, mse, singleEstimates)
	if(showEmpiricalMse)
	{
		retVal <- results[indices,c("method", "sampleSize", "wnrv", "relativeErrors", "mse", "singleEstimates")]
		colnames(retVal) <- c("Method", "$n$", "WNRV", "RE", "MSE", "$\\widehat{\\ell}$")
	}
	else
	{
		retVal <- results[indices,c("method", "sampleSize", "wnrv", "relativeErrors", "singleEstimates")]
		colnames(retVal) <- c("Method", "$n$", "WNRV", "RE", "$\\widehat{\\ell}$")
	}
	retVal$Method <- c("approxZeroVarianceWOR" = "WOR", "approxZeroVariance" = "IS")[retVal$Method]
	return(retVal)
}
indices <- which(scenarios$graph == "dodecahedron" & scenarios$probability == 0.99)
print(xtable(produceTable(indices, TRUE), digits = c(0, 0, 0, 2, 2, 2, 6), display = c("s", "s", "d", "e", "e", "e", "e"), label = "table:networkReliability_dodec", caption = "Simulation results for the dodecahedron graph with $p = 0.99$."), include.rownames=FALSE, math.style.exponents = TRUE, sanitize.text.function = identity)

indices <- which(scenarios$graph == "grid9Augmented" & scenarios$probability == 0.9999)
print(xtable(produceTable(indices, FALSE), digits = c(0, 0, 0, 2, 2, 6), display = c("s", "s", "d", "e", "e", "e"), label = "table:networkReliability_grid9Augmented", caption = "Simulation results for the augmented $9\\times 9$ grid graph with $p = 0.9999$."), include.rownames=FALSE, math.style.exponents = TRUE, sanitize.text.function = identity)

indices <- which(scenarios$graph == "dodecParallel" & scenarios$probability == 0.9999)
print(xtable(produceTable(indices, FALSE), digits = c(0, 0, 0, 2, 2, 6), display = c("s", "s", "d", "e", "e", "e"), label = "table:networkReliability_dodecParallel", caption = "Simulation results for the three dodecahedron graphs in parallel with $p = 0.9999$."), include.rownames=FALSE, math.style.exponents = TRUE, sanitize.text.function = identity)

indices <- which(scenarios$graph == "dodecSeries" & scenarios$probability == 0.9999)
print(xtable(produceTable(indices, FALSE), digits = c(0, 0, 0, 2, 2, 6), display = c("s", "s", "d", "e", "e", "e"), label = "table:networkReliability_dodecSeries", caption = "Simulation results for the three dodecahedron graphs in series with $p = 0.9999$."), include.rownames=FALSE, math.style.exponents = TRUE, sanitize.text.function = identity)
