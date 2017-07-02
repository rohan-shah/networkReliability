library(boot)
library(networkReliability)
library(ggplot2)
library(Hmisc)
library(tikzDevice)

load("./results/approxZeroVarianceWORWithVariance-20-dodecahedron-0.9999.RData")

empiricalResults <- unlist(lapply(results, function(x) as.numeric(x@estimate)))
htResults <- unlist(lapply(results, function(x) as.numeric(x@varianceEstimate)))

set.seed(1)
nSamples <- 20000
sampleSize <- 1000
bootstrapResults <- lapply(1:nSamples, function(x)
{
        var(empiricalResults[sample(1:length(empiricalResults), sampleSize)])
})

tikz("./htEstimates.tex", standAlone = TRUE)
breaks <- seq(1e-32, 3.5e-29, length.out = 40)
labelBreaks <- c(breaks[seq(1, length(breaks), 10)], tail(breaks, 1))
labelBreaksLatex <- paste0("$", latexSN(signif(labelBreaks, 2)), "$")
print(ggplot(data.frame(value = htResults), aes(value)) + stat_bin(breaks = breaks, pad = FALSE) + scale_x_continuous(breaks = labelBreaks, labels = labelBreaksLatex) + xlab("Horvitz-Thompson Variance") + ylab("Count") + ggtitle("") + theme(axis.text = element_text(size = rel(1.35)), axis.title = element_text(size = rel(2)), plot.margin = margin(0, 25, 5, 5)) + geom_vline(xintercept = var(empiricalResults), size = 2) + geom_vline(xintercept = mean(htResults), size = 2, color = "red"))
dev.off()

tikz("./bootstrapEmpirical.tex", standAlone = TRUE)
print(ggplot(data.frame(value = log10(unlist(bootstrapResults))), aes(value)) + stat_bin() + xlab("$\\log_{10}\\left(\\mathrm{Bootstrapped Variance}\\right)$") + ylab("Count") + ggtitle("") + theme(axis.text = element_text(size = rel(1.35)), axis.title = element_text(size = rel(2)), plot.margin = margin(0, 25, 5, 5)) + geom_vline(xintercept = log10(var(empiricalResults)), size = 2) + geom_vline(xintercept = log10(mean(htResults)), size = 2, color = "red"))
dev.off()
