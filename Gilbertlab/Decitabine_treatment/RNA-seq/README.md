# RNA-seq data 

## Analysis workflow 
**Differential RNA Expression:** I've used `salmon` for alignment and `DESeq2` for ANNOVA like analysis (see [DESeq2 for time-series-experiments](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments.)). Mainly, the linear model designed as `~condition + time + condition:time` to include both time and treatment together in the model. 
- See [HL60 cell line results](https://github.com/goodarzilab/Abe/blob/master/Gilbertlab/Decitabine_treatment/RNA-seq/hl60-exp/README.md)

**Differential RNA Stability:** Separately, I've used `STAR` for alignment and then, I've used `featureCount` to evalute intronic and exonic counts (similar to [CRIES](https://github.com/csglab/CRIES)). Then, [REMBRANDTS](https://github.com/csglab/REMBRANDTS) is the tool to estimate unbias transcript stability. Finally, I use `limma` for differential analysis. 
- See [HL60 cell line results]https://github.com/goodarzilab/Abe/blob/master/Gilbertlab/Decitabine_treatment/RNA-seq/hl60-stbl/README.md)
