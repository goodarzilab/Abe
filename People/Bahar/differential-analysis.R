args <- commandArgs(trailingOnly = TRUE)

countsDIR <- args[1]
colDataDIR <- args[2]
name <- args[3]
# PDIR <- args[1]
# refCOND <- [4]
# JOBS <- args[5]
# setwd(PDIR)
dir.create(paste0('deseq/',name))


suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(patchwork))
# suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(SingleCellExperiment))


sce_prep <- function (counts,colData){
    ### Load data as `SingleCellExperiment` object 
    sce <- SingleCellExperiment(list(counts=counts),colData=colData)
    
    # ## low count filter - at least 10 cells with count of 5 or more
    # keep <- rowSums(counts(sce) >= 1) >= 10 
    # sce <- sce[keep,] 
    
    assay(sce) <- as.matrix(assay(sce))

    return (sce)
}

plotDisp <- function(dds){
    p1 = plotDispEsts(dds)
    keep <- rowSums(counts(dds) >= 10) >= 25
    dds2 <- estimateDispersionsFit(dds[keepForDispTrend,])
    p2 = plotDispEsts(dds2, ylim=c(1e-3,1))
}

plot_PCA <- function(vsd, gr, title,legend="none"){
    colData <- colData(vsd)
    z <- plotPCA(vsd,intgroup=gr, returnData=TRUE)
    percentVar <- round(100 * attr(z, "percentVar"))
    pca <- ggplot(z, aes(PC1, PC2)) + 
            geom_point(aes(size = 0.01, colour=group,fill=group), alpha = 0.1,shape=18) +
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
            ylab(paste0("PC2: ",percentVar[2],"% variance")) +
            guides (size = FALSE) +
            ggtitle (title)+ 
            theme(legend.position=legend)
    return (pca)
}

correct_batch <- function (dds,gr,out){
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    p0 <- plot_PCA(vsd, gr, 'Before removeBatchEffect')
    mat <- assay(vsd)
    
    mat <- limma::removeBatchEffect(mat, vsd$batch) 
    assay(vsd) <- mat
    p1 <- plot_PCA(vsd, gr, 'After removeBatchEffect',legend="right")
    counts_batch_corrected <- assay(vsd)
    
    if (out == 'plot'){return (p0 + p1)}
    if (out == 'vsd') {return (vsd)}
    if (out == 'cbc') {return (counts_batch_corrected)}
} 

get_res<- function(fit,coef){
    fc <- fit[,coef]$coef [,coef] %>% data.frame 
    fc <- fc %>% dplyr::rename(log2FoldChange=colnames(fc))
    pv <- fit[,coef]$p.value [,coef] %>% data.frame 
    pv <- pv %>% dplyr::rename(pvalue=colnames(pv))
    res <- cbind(fc, pv)
    return(res)
}


###################################################################################################
message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
message ('Part 1. ')
message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
###################################################################################################

counts <- fread(countsDIR,sep='\t') %>% 
    data.frame 
counts <- counts[!duplicated(counts$V1),] %>% remove_rownames %>% 
    column_to_rownames('V1') %>% 
    # as.matrix %>% 
    mutate_all(function(x) as.numeric(as.character(x)))

colData <- fread(colDataDIR,sep='\t') %>% 
    data.frame %>% 
    column_to_rownames('V1')

rownames(colData) <- colnames(counts) 

colData$rep = factor(colData$rep)

colData$cond <- factor (colData$cond)
colData$cond <- relevel(colData$cond, ref="Balb")

conds <- colData$cond 
reps  <- colData$rep

print ('counts loaded!')
print (counts %>% dim)
print ('colData loaded!')
print (colData %>% dim)

###################################################################################################
message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
message ('Part 2. ')
message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
###################################################################################################
sc <- sce_prep(counts,colData)
de <- convertTo(sc, type="edgeR")

message ('object created!')
message (class(de))

design <- model.matrix(~0 + conds + reps )
colnames(design) <- gsub("conds", "", colnames(design))

contr.matrix <- makeContrasts(
    NSGvsBalb= NSG-Balb, 
    NujvsBalb= Nuj-Balb, 
    RagvsBalb= Rag-Balb, 
    NSGvsRag = NSG-Rag, 
    NujvsRag = Nuj-Rag, 
    levels = colnames(design)
)

message ('Model:')
message (contr.matrix)

# Log tranform & Quantile normalize
de2 <- normalizeBetweenArrays(log(de$counts + 0.001))

vfit2 <- lmFit(de2, design)
vfit2 <- contrasts.fit(vfit2, contrasts=contr.matrix)
efit2 <- eBayes(vfit2)
# plotSA(efit, main="Final model: Mean-variance trend")

res_NSGvsBalb2 <- get_res(efit2,'NSGvsBalb')
write.table(
    res_NSGvsBalb2,
    paste0('deseq/',name,'/delta_exp_NSG_vs_Balb.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)

res_NujvsBalb2 <- get_res(efit2,'NujvsBalb')
write.table(
    res_NujvsBalb2,
    paste0('deseq/',name,'/delta_exp_Nuj_vs_Balb.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)

res_RagvsBalb2 <- get_res(efit2,'RagvsBalb')
write.table(
    res_RagvsBalb2,
    paste0('deseq/',name,'/delta_exp_Rag_vs_Balb.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)

res_NSGvsRag2  <- get_res(efit2,'NSGvsRag')
write.table(
    res_NSGvsRag2,
    paste0('deseq/',name,'/delta_exp_NSG_vs_Rag.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)

res_NujvsRag2  <- get_res(efit2,'NujvsRag')
write.table(
    res_NujvsRag2,
    paste0('deseq/',name,'/delta_exp_Nuj_vs_Rag.txt'),
    sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE
)


    
    
###################################################################################################
message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
message ('Part 3. ')
message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
###################################################################################################
plot_Volcano <- function(res, lfc.cutoff  = 1,pval.cutoff = 0.05, title=''){#, x_min=-20,x_max=20, y_max=FALSE){
    res$sig <- as.factor(res$pvalue < pval.cutoff & abs(res$log2FoldChange) > lfc.cutoff)
    relevel(res$sig, ref=TRUE)
    res <- res %>% mutate(name=rownames(res))
    
    data = subset(res[order(res$pvalue),], sig == TRUE)[1:7,]
    
    # if (y_max == FALSE){
    -log10(min(res$pvalue))
    # }
    
    vol = res %>% ggplot(
        aes(x=log2FoldChange, y=-log10(pvalue), colour=sig, fill=sig)) +
        geom_point(aes(color = sig),alpha = 1/10) +
            geom_hline(yintercept=-log10(pval.cutoff), linetype="dashed", alpha = 4/10) +
            geom_vline(xintercept=lfc.cutoff, linetype="dashed", alpha = 4/10) +
            geom_vline(xintercept=(-1)*lfc.cutoff, linetype="dashed", alpha = 4/10) +
            scale_color_manual(values = c("grey", "red")) +
            theme_bw() + 
            theme(
                legend.position="none",plot.title = element_text(hjust = 0.5)
            ) +
            ggtitle (title) + 
            ggrepel::geom_text_repel(
                data = data,
                aes(label = rownames(data)),
                size = 3,
                box.padding = unit(0.35, "lines"),
                point.padding = unit(0.3, "lines")
            ) 
            # + 
            # xlim(c(x_min,x_max)) +
            # if (y_max != FALSE){ylim(c(0,y_max))}

     return (vol)
} 


# y_max=250
# x=5

p1 <- res_NSGvsBalb2 %>% plot_Volcano(
    lfc.cutoff  = 0.1,
    title='NSGvsBalb'
) /
res_NujvsBalb2 %>% plot_Volcano(
    lfc.cutoff  = 0.1,
    title='NujvsBalb'
) /
res_RagvsBalb2 %>% plot_Volcano(
    lfc.cutoff  = 0.1,
    title='RagvsBalb'
) 


p2 <- res_NSGvsRag2  %>% plot_Volcano(
    lfc.cutoff  = 0.1,
    title='NSGvsRug'
) /
res_NujvsRag2  %>% plot_Volcano(
    lfc.cutoff  = 0.1,
    title='NujvsRag'
) 
p2

ggsave(paste0('deseq/',name,'/volcano_1.pdf'),p1)
ggsave(paste0('deseq/',name,'/volcano_2.pdf'),p2)

message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
message ('@@@@@@@@ DONE! :-)@@@@@@@@@@@')
message ('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
message (date())
# sessionInfo()
