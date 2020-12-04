library (exomePeak) 
run_exomepeak <- function(bams, enz){    
    TREATED <- bams[-grep('NT', bams)]
    nonTREATED <- bams[grep('NT', bams)]
    # Design: 
    IP_BAM= nonTREATED[grep('meRIP',nonTREATED)]
    INPUT_BAM=nonTREATED[grep('IN',nonTREATED)]
    # get only bams TREATED with input enz (for experiments with multiple enz treatments)
    T <- TREATED[grep(enz,TREATED)]
    TREATED_IP_BAM <- T[grep('meRIP',T)]
    TREATED_INPUT_BAM <- T[grep('IN',T)]
    # comparison 
    
    res <- exomepeak(GENOME="hg38", 
                     IP_BAM=IP_BAM, 
                     INPUT_BAM=INPUT_BAM, 
                     TREATED_IP_BAM=TREATED_IP_BAM, 
                     TREATED_INPUT_BAM=TREATED_INPUT_BAM)
    return (res)
}