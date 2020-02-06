library(tidyverse)
library(argparse)
library(GenomicFeatures)
library(exomePeak)

parser <- ArgumentParser()
parser$add_argument('-i', '--input_dir', nargs=1, type='character', help='Absolute path to the input directories')
parser$add_argument('-s', '--species', nargs=1, type='character', help='Species hg19, hg38 (default), or mm10', default='hg38')
parser$add_argument('-e', '--enzymes', nargs=1, type='character', help='list of enzymes separated by ,')
parser$add_argument('-m', '--mode', nargs=1, type='character', help='select run mode (control, treated)')
args <- parser$parse_args()
indir <- args$input_dir
setwd(indir)

######################################## functions ########################################
if (args$species == 'hg38') {
  txdb <- makeTxDbFromGFF('/rumi/shams/genomes/hg38/Homo_sapiens.GRCh38.87.gtf',
                          organism='Homo sapiens')
  tag <- '.bam$'
#} else if (args$species == 'mm10') {
  #
} else if (args$species == 'HIV') {
  txdb <- makeTxDbFromGFF('./HIV/HIV.gtf')
  tag <- '.HIV.bam$'
}
setwd("./bam")
print (txdb)

if (args$mode == 'treated'){
  run_exomepeak <- function (enz){
    res <- exomepeak(
      TXDB=txdb,
      IP_BAM=grep('_meRIP_', grep(paste('NT',tag,sep=''), list.files(), value=1) , value=1),
      INPUT_BAM=grep('_IN_', grep(paste('NT',tag,sep=''), list.files(), value=1) , value=1),
      TREATED_IP_BAM=grep('_meRIP_',grep(paste(enz,tag,sep=''),list.files(), value=1), value=1),
      TREATED_INPUT_BAM=grep('_IN_',grep(paste(enz,tag,sep=''),list.files(), value=1), value=1),
      OUTPUT_DIR=paste('../exomepeak/', args$species, sep=''),
      EXPERIMENT_NAME=enz
    )
    saveRDS(res, paste('../exomepeak',args$species, enz,'results.rds', sep='/'))
    return (res)
}} else if (args$mode == 'control') {
  run_exomepeak <- function (){
    res <- exomepeak(
      TXDB=txdb,
      IP_BAM=grep('_meRIP_', grep(paste('NT',tag,sep=''), list.files(), value=1) , value=1),
      INPUT_BAM=grep('_IN_', grep(paste('NT',tag,sep=''), list.files(), value=1) , value=1),
      OUTPUT_DIR=paste('../exomepeak/', args$species, sep=''),
      EXPERIMENT_NAME='control'
    )
    saveRDS(res, paste('../exomepeak',args$species, 'control', 'results.rds', sep='/'))
    return (res)
}}
######################################## run script ########################################
enzymes <- unlist(str_split(args$enzymes, ','))
if (args$mode == 'treated'){
  for (enz in enzymes){
    run_exomepeak(enz)
  }
}
if (args$mode == 'control'){
  run_exomepeak()
}
