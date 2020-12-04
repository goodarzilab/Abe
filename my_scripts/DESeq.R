library (GenomicFeatures)
library (tximport)
library (tidyverse)
library (ggplot2)
library (ggrepel)
library (DESeq2)
library (patchwork)
library (BiocParallel)
register(MulticoreParam(4))

"""
part 0: useful functions 

"""
Box_plot <- function (dds, name_it){
    pdf(paste(name_it,'pdf',sep='.') )
    b = boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
    dev.off() 

    jpeg(paste(name_it,'jpg',sep='.') )
    b = boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
    dev.off()
}

ggplot_Save <- function (p, name_it){
    ggsave(paste(name_it,'png',sep='.'), plot = p, device = 'png', dpi = 300)
    ggsave(paste(name_it,'pdf',sep='.'), plot = p, device = 'pdf', dpi = 300)
}

plot_PCA <- function(dds, title=''){
    vsd <- varianceStabilizingTransformation(dds)
    z <- plotPCA(vsd,intgroup=c('Myc','Met'), returnData=TRUE)
    percentVar <- round(100 * attr(z, "percentVar"))

    pca <- ggplot(z, aes(PC1, PC2)) +
            geom_point(aes(size = 2,  shape=group), alpha = 4/10) +
            geom_text_repel(aes(label = row.names(colData)),size = 3.5) +
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
            ylab(paste0("PC2: ",percentVar[2],"% variance")) +
            guides (size = FALSE) +
            ggtitle (title)+ 
            theme(legend.position="none")
    return (pca)
}

get_result <- function (dds, design){
    res <- results(dds, contrast=design, parallel=TRUE) %>% data.frame %>% rownames_to_column('gene_id')
    gn2ens = getBM(
        attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
        filters = 'ensembl_gene_id_version', 
        values = res$gene_id,
        mart = ensembl
    ) %>%  data.frame() %>% column_to_rownames('ensembl_gene_id_version')
    
    res <- res %>% column_to_rownames('gene_id')
    return (res)
}

plot_Volcano <- function(res, lfc.cutoff  = 1,pval.cutoff = 0.05, title=''){
    res$sig <- as.factor(res$pvalue < pval.cutoff & abs(res$log2FoldChange) > lfc.cutoff)
    relevel(res$sig, ref=TRUE)
    ## Volcano plot
    vol = res %>% ggplot(
        aes(x=log2FoldChange, y=-log10(pvalue), colour=sig, fill=sig)) +
        geom_point(aes(color = sig),alpha = 1/10) +
    #         xlim(c(-20,20)) +
    #         ylim(c(0,11)) +
            geom_hline(yintercept=-log10(pval.cutoff), linetype="dashed", alpha = 4/10) +
            geom_vline(xintercept=lfc.cutoff, linetype="dashed", alpha = 4/10) +
            geom_vline(xintercept=(-1)*lfc.cutoff, linetype="dashed", alpha = 4/10) +
            scale_color_manual(values = c("grey", "red")) +
            theme_bw() + 
            theme(legend.position="none") +
            ggtitle (title) + 
            geom_text_repel(
                data = subset(res[order(res$pvalue),], sig == TRUE)[1:5,],
                aes(label = name),
                size = 3,
                box.padding = unit(0.35, "lines"),
                point.padding = unit(0.3, "lines")
            )
     return (vol)
}

'''
part 1: read annotations
'''
# https://www.biostars.org/p/196367/
# This is the final solution!!
gtf <- rtracklayer::import('/rumi/shams/abe/genomes/hg38/gencode.v28/gencode.v28.annotation.gtf')
gene2name <- gtf[gtf$type == "gene"] %>% data.frame %>% column_to_rownames('gene_id') %>% dplyr::select('gene_name')

txdb  = makeTxDbFromGFF(
    '/rumi/shams/abe/genomes/hg38/gencode.v28/gencode.v28.annotation.gtf', 
    organism='Homo sapiens')

# tx2gene objects 
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")


'''
part 2: read  data
'''
files <- list.files(path='./quants', pattern="quant.sf",full.names = TRUE, recursive=T)
names(files) <- gsub("./quants/(\\S+)/quant.sf","\\1",files)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut=T)

txi.gene <- summarizeToGene(txi, tx2gene, ignoreAfterBar= TRUE)

# D is lowMyc
# Dpo is highMyc
# MF is highMet
# M10 is lowMet
Myc  = c(rep('lowMyc', 6), rep('highMyc', 6) )
Met  = rep(c(rep('highMet',3), rep('lowMet',3) ),2)
reps = rep(c('1','2','3'), 4)

colData <- data.frame(Myc = Myc, Met = Met,sample_id=paste(Myc,Met,reps, sep='_'),row.names=colnames(txi$abundance))

'''
part 3: Principal component analysis
'''
dds0 <- DESeqDataSetFromTximport(txi, colData, ~Myc + Met )
dds_raw <- DESeq(dds0, parallel=TRUE)
# boxplot
Box_plot(dds_raw, 'boxplot_raw')
# filtering 
keep <- rowSums(counts(dds0) > 10) == ncol(dds0)
dds1 <- dds0[keep,]
dds_filter <- DESeq(dds1, parallel=TRUE)
# boxplot
Box_plot(dds_filter, 'boxplot_filter')
# PCA 
pca = plot_PCA(dds_filter)
ggplot_Save(pca, 'PCA')

'''
part 4: run the DESeq tests
'''
dds <- DESeqDataSetFromTximport(txi.gene, colData, ~Myc + Met)
dds <- DESeq(dds, parallel=TRUE)
ncu <- counts(dds, normalized=TRUE)
write.table(ncu,'DE2norm.txt', sep="\t", quote=FALSE, col.names=paste(names(files), colData$sample_id, sep=':') )

res_Met <- results(dds, contrast=list('Met_lowMet_vs_highMet'), parallel=TRUE) %>% data.frame %>% add_column(name=gn2name[rownames(dds),])
write.table(res_Met,'Met_DE2res.txt', sep="\t", quote=FALSE, col.names=NA)
p = plot_Volcano (res_Met, title='lowMet_vs_highMet')
ggplot_Save(p, 'Met_Volcano')

res_Myc <- results(dds, contrast=list('Myc_lowMyc_vs_highMyc'), parallel=TRUE) %>% data.frame %>% add_column(name=gn2name[rownames(dds),])

write.table(res_Myc,'Myc_DE2res.txt', sep="\t", quote=FALSE, col.names=NA)
p = plot_Volcano (res_Myc, title='lowMyc_vs_highMyc')
ggplot_Save(p, 'Myc_Volcano')

