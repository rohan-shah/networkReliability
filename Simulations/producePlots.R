load("./summarised.RData")
source("./generateScenarios.R")
library(ggplot2)
library(scales)
library(dplyr)
library(scales)
theme_set(theme_bw() + theme(plot.margin = unit(c(5.5, 12, 5.5, 5.5), "points")))
size <- 1.5

data <- cbind(scenarios, wnrv = wnrv, relativeError = relativeErrors, empiricalVariance = empiricalVariances, workNormalizedVariance = workNormalizedVariance)
rewrites <- c("approxZeroVariance" = "IS", "approxZeroVarianceWOR" = "WOR", "approxZeroVarianceWORMerge" = "WOR-Merge", "fearnhead" = "Fearnhead")
data <- subset(data, method %in% names(rewrites))
data$method <- rewrites[data$method]

scale_y <- scale_y_log10(label = trans_format("log10", math_format(10^.x)))
formatTheme <- theme(axis.text.x = element_text(size = 22), axis.title.x = element_text(size = 24), axis.text.y = element_text(size = 22), axis.title.y = element_text(size = 20), legend.text = element_text(size = 18), legend.title = element_text(size = 18), legend.key.width = unit(4, "line"))
custom_linetype <- scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash"))

plot <- data %>% filter(graph == "dodecahedron" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("WNV") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./dodecahedronWNV0_99.pdf")
print(plot)
dev.off()
#Plot RE, for the dodecahedron with p = 0.99
plot <- data %>% filter(graph == "dodecahedron" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("RE") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./dodecahedronRE0_99.pdf")
print(plot)
dev.off()


#Plot workNormalizedVariance, for the augmented grid graph with p = 0.99
plot <- data %>% filter(graph == "grid9Augmented" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("WNV") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./grid9AugmentedWNV0_99.pdf")
print(plot)
dev.off()
#Plot RE, for the augmented grid graph with p = 0.99
plot <- data %>% filter(graph == "grid9Augmented" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("RE") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./grid9AugmentedRE0_99.pdf")
print(plot)
dev.off()

#Plot workNormalizedVariance, for the grid graph with p = 0.99
plot <- data %>% filter(graph == "grid9" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("WNV") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./grid9WNV0_99.pdf")
print(plot)
dev.off()
#Plot RE, for the grid graph with p = 0.99
plot <- data %>% filter(graph == "grid9" & probability == 0.99 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("RE") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./grid9RE0_99.pdf")
print(plot)
dev.off()

#Plot workNormalizedVariance, for the dodecahedron series graph with p = 0.9999
plot <- data %>% filter(graph == "dodecSeries" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("WNV") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./dodecSeriesWNV0_9999.pdf")
print(plot)
dev.off()
#Plot RE, for the dodecahedron series graph with p = 0.9999
plot <- data %>% filter(graph == "dodecSeries" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("RE") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./dodecSeriesRE0_9999.pdf")
print(plot)
dev.off()
#Plot variance, for the dodecahedron series graph with p = 0.9999
plot <- data %>% filter(graph == "dodecSeries" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, empiricalVariance, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("Var") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./dodecSeriesVar0_9999.pdf")
print(plot)
dev.off()

#Plot workNormalizedVariance, for the dodecahedron parallel graph with p = 0.9999
plot <- data %>% filter(graph == "dodecParallel" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, workNormalizedVariance, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("WNV") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./dodecParallelWNV0_9999.pdf")
print(plot)
dev.off()
#Plot RE, for the dodecahedron parallel graph with p = 0.9999
plot <- data %>% filter(graph == "dodecParallel" & probability == 0.9999 & method != "approxZeroVarianceWORWithVariance" & method != "residualResampling") %>% ggplot(mapping = aes(sampleSize, relativeError, colour = method, linetype = method)) + scale_y + geom_line(size = size) + scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = c(10, 100, 1000, 10000)) + theme(legend.position = c(0.3, 0.2)) + formatTheme + xlab("Sample size") + ylab("RE") + guides(colour = guide_legend(title = "Method"), linetype = guide_legend(title = "Method")) + custom_linetype
pdf("./dodecParallelRE0_9999.pdf")
print(plot)
dev.off()
