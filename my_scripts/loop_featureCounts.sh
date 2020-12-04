#!/usr/bin/env bash
#loop featureCounts intron and exon counts
bam_folder=$1
stbl_folder=$2

set -o errexit
set -o nounset
set -o xtrace

GTF_index='/rumi/shams/genomes/hg38/hg38_ensemble_'

main() {
  local f
  for f in ${bam_folder}/*.bam; do
    base=`basename $f`
    sample=${base/.bam/};
    echo -e '----------------------- ' $sample  ' -----------------------'
    echo featureCounts -M -T 18 -t intron -g gene_id -a ${GTF_index}introns.gtf -o ${stbl_folder}/counts/${sample}_introns.txt ${f}
    echo featureCounts -M -T 18 -t exon -g gene_id -a ${GTF_index}consExons.gtf -o ${stbl_folder}/counts/${sample}_exons.txt ${f}
  done
}

mkdir -p $stbl_folder
mkdir -p ${stbl_folder}/counts

main
