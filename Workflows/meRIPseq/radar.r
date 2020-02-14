library("RADAR")

parser <- ArgumentParser()
parser$add_argument('-i', '--input_dir', nargs=1, type='character', help='Absolute path to the input directories')
parser$add_argument('-g', '--gtffile', nargs=1, type='character', help='Absolute path to the gtf file')
parser$add_argument('-c', '--comparisons', nargs=1, type='character', help='List of comparisons between conditions: BvsA,CvsB... where A, B... are baselines, respectively')
parser$add_argument('-ip', '--modification', nargs=1, type='character', help='modification')
parser$add_argument('-s', '--samples', nargs=1, type='character', help='Sample prefix in file names prior to treatment and INPUT/IP')
parser$add_argument('-t', '--threads', nargs=1, help='Number of the threads. 16 by default', type='integer', default=16)

args <- parser$parse_args()
indir <- args$input_dir
gtffile = args$gtffile

meta <- read.table(args$meta, sep='\t', header=FALSE)
comp = as.data.frame(t(data.frame(strsplit(strsplit(args$comparisons, ',')[[1]], 'vs'), row.names=c('condition_1','condition_2'))))

setwd(args$input_dir)
setwd('./bam')



run_radar <- function(enz, gtffile,species, cutoff = 0.1, Beta_cutoff = 0.5,threads = 18){
    outputDir = paste("radar", species,enz, sep='/')
    radar <- countReads(
        samplenames = unlist(lapply (c('s23','s24'), paste, paste(species,c('NT',enz), sep='.'),sep='.')),
        gtf = gtffile,
        bamFolder = "bam",
        modification = 'm6A',
        strandToKeep = "opposite",
        outputDir = outputDir,
        threads = threads,
        saveOutput = TRUE
    )
    radar <- normalizeLibrary(radar) #, boxPlot = FALSE)
    radar <- adjustExprLevel(radar)
    variable(radar) <- data.frame( Group =data.frame( Group = rep(c("input","meRIP"),2)) )
    radar <- filterBins(radar,minCountsCutOff = 15)
    radar <- diffIP_parallel(radar, thread = threads)
    top_bins <- extractIP(radar,filtered = T)[order(rowMeans( extractIP(radar,filtered = T) ),decreasing = T)[1:1000],]
    radar <- reportResult(radar, cutoff = cutoff, Beta_cutoff = Beta_cutoff, threads=threads)
    result <- results(radar)
    saveRDS(radar, file = paste(outputDir,"radar.rds",sep='/'))
}
