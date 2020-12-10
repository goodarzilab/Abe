#!/usr/bin/env bash
#loop featureCounts intron and exon counts

set -o errexit
set -o nounset
set -o xtrace

bam_folder='/rumi/shams/abe/Gilbertlab/Decitabine_treatment/RNA-seq/hl60-bam'
stbl_folder='/rumi/shams/abe/Gilbertlab/Decitabine_treatment/RNA-seq/hl60-stbl'
GTF_index='/rumi/shams/genomes/hg38/hg38_ensemble_'

main() {
  local f
  for f in ${bam_folder}/*.bam; do
    base=`basename $f`
    sample=${base/.bam/};
    echo -e '----------------------- ' $sample  ' -----------------------'
    featureCounts -M -T 18 -t intron -g gene_id -a ${GTF_index}introns.gtf -o ${stbl_folder}/counts/${sample}_introns.txt ${f}
    featureCounts -M -T 18 -t exon -g gene_id -a ${GTF_index}consExons.gtf -o ${stbl_folder}/counts/${sample}_exons.txt ${f}
  done
}

mkdir -p $stbl_folder
mkdir -p ${stbl_folder}/counts

main
