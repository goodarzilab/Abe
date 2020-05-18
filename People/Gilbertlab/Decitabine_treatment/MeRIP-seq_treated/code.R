library(RADAR)
setwd('~/People/Gilbertlab/Decitabine_treatment/MeRIP-seq_treated/')
# load bam files 
radar <- countReads(samplenames=c('T1','T2','U1','U2'), gtf = "/rumi/shams/genomes/hg38/hg38_genes.gtf",bamFolder='./bam/',modification = "m6A", strandToKeep = "opposite",outputDir='radar', threads = 18)

saveRDS(radar, 'radar/raw.rds')

pdf('radar/boxPlot.pdf')
radar <- normalizeLibrary(radar)
dev.off()
png('radar/boxPlot.png')
radar <- normalizeLibrary(radar)
dev.off()

radar <- adjustExprLevel(radar)

variable(radar) <- data.frame( Group = c(rep("T",2),rep("U",2)) )

radar <- filterBins(radar,minCountsCutOff = 15)
radar <- diffIP_parallel(radar, thread = 8)
top_bins <- extractIP(radar,filtered = T)[order(rowMeans( extractIP(radar,filtered = T) ),decreasing = T)[1:1000],]

png('radar/PCAPlot.png')
plotPCAfromMatrix(top_bins,group = unlist(variable(radar)) )
dev.off()
pdf('radar/PCAPlot.pdf')
plotPCAfromMatrix(top_bins,group = unlist(variable(radar)) )
plotPCAfromMatrix(top_bins,group = unlist(variable(radar)) )
dev.off()

radar <- reportResult(radar, cutoff = 0.1, Beta_cutoff = 0.5)
saveRDS(radar, 'radar/final.rds')

result <- results(radar)

pdf('radar/HeatmapPlot.pdf')
plotHeatMap(radar)
dev.off()
png('radar/HeatmapPlot.png')
plotHeatMap(radar)
dev.off()

png('radar/peakDistribution.png')
peakDistribution(radar)
dev.off()

res = results(radar)

w <- wilcox.test( res$logFC, mu=0, alternative = "greater")
t <- t.test(      res$logFC, mu=0, alternative = "greater")
h = ggplot(res, aes(x=logFC)) + 
    ggtitle('Decitabine_treatment',
sprintf("wilcox.test (p.value): %.5f \t t.test (p.value):%.5f \t [mu=0,alter=greater]", w$p.value, t$p.value)) +geom_histogram(binwidth=0.1)
ggsave("radar/peak_Histograms.png", plot = h, device = 'png', dpi = 300)
ggsave("radar/peak_Histograms.pdf", plot = h, device = 'pdf', dpi = 300)


v = res %>% ggplot(aes(x=logFC, y=-log10(p_value),colour=sig, fill=sig)) +geom_point(aes(color = sig),alpha = 1/10) +ggtitle('Decitabine_treatment') +
geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha = 4/10) +
geom_vline(xintercept=-3, linetype="dashed", alpha = 4/10) +
geom_vline(xintercept= 3, linetype="dashed", alpha = 4/10) +
scale_color_manual(values = c("grey", "red")) +
theme_bw() + theme(legend.position="none") +
geom_text_repel(
data = subset(res[order(res$p_value),], sig == TRUE)[1:10,],aes(label = name),size = 5,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))
ggsave("radar/peak_Volcanos.png", plot = v, device = 'png', dpi = 300)
ggsave("radar/peak_Volcanos.pdf", plot = v, device = 'pdf', dpi = 300)

radar <- PrepCoveragePlot(radar)

png('radar/top_peak.png')
plotGeneCov(radar,geneName = "ZNF93", center = median, libraryType = "opposite")
dev.off()

