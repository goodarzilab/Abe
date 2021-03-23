mkdir -p bam
mkdir -p qc_star

STAR --genomeLoad LoadAndExit --genomeDir ~/genomes/hg38/gencode.v34/star_index

for fq in fastq/*R1*; do
    fq=`basename $fq`
    echo '------------' $fq '-----------'
    out=${fq/_S*/}
    STAR \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --runThreadN 18 \
    --genomeDir ~/genomes/hg38/gencode.v34/star_index \
    --readFilesIn fastq/$fq \
    --outFileNamePrefix bam/$out
done

STAR --genomeLoad Remove --genomeDir ~/genomes/hg38/gencode.v34/star_index

rm -vr _STARtmp/ Aligned.out.sam Log.out Log.progress.out

