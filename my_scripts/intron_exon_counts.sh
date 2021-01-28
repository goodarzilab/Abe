#loop featureCounts intron and exon counts

PDIR=$1
bamDIR=$2
stblDIR=$3
JOBS=$4

GTF_index='/rumi/shams/genomes/hg38/hg38_ensemble_'

mkdir -p $stblDIR
mkdir -p ${stblDIR}/counts
mkdir -p ${stblDIR}/logs

for f in ${bam_folder}/*.bam; do
    base=`basename $f`
    sample=${base/.bam/};
    echo -e '----------------------- ' $sample  ' -----------------------'
    featureCounts -M -T $JOBS -t intron -g gene_id -a ${GTF_index}introns.gtf -o ${stblDIR}/counts/${sample}_introns.txt ${f} &> ${stblDIR}/logs/${sample}_introns.log
    featureCounts -M -T $JOBS -t exon -g gene_id -a ${GTF_index}consExons.gtf -o ${stblDIR}/counts/${sample}_exons.txt ${f} &> ${stblDIR}/logs/${sample}_exons.log
done
