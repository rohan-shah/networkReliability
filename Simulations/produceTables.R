source("./generateScenarios.R")
library(xtable)
load("./summarised.RData")

indices <- which(scenarios$method %in% "approxZeroVarianceWORWithVariance")
data <- cbind(scenarios, varianceEstimates, empiricalVariances)[indices,c("graph", "sampleSize", "probability", "varianceEstimates", "empiricalVariances")]
print(xtable(data, digits = c(0, 0, 0, 4, 3, 3), display = c("s", "s", "d", "f", "e", "e"), label = "table:network_reliability_variances", caption = "Simulation results for the Horvitz-Thompson estimate of the variance."), include.rownames=FALSE, math.style.exponents = TRUE, sanitize.text.function = identity)
