# ############################### Trimming ############################################# #
# for f in *.fastq.gz; do
#     out=${f/R1_001.fastq.gz/trim.fastq.gz};     
#     echo cutadapt -u 3 -j 12 -q 15 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o $out $f;
# done

# ############################### Alignment ############################################ #
# for f in fastq/*trim.fastq.gz;
# 	do     out=${f/_S[0-9]*/_};
#	out=${out/fastq/bam};
#	STAR --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN 16 --sjdbGTFfile /rumi/shams/genomes/hg38/Homo_sapiens.GRCh38.87.gtf --genomeDir /rumi/shams/genomes/hg38 --readFilesIn $f --outFileNamePrefix $out --outReadsUnmapped Fastx; 
# done

# # ############################### Alignment: HIV genome ############################## #
# for f in fastq_unmapped/*; 
# 	do o=${f/_Unmapped.out.mate1/.HIV.bam}; 
# 	o=${o/fastq_unmapped/bam}; 
# 	bowtie2 --sensitive -N 1 -x HIV/HIV -U $f | samtools sort -o $o; 
# done

# # ############################### counts ############################################# #
# for f in bam/*_*.bam; do
#     out=${f/.bam/.fc};
#     out=${out/bam/fc};
#     featureCounts -T 12 -P -B -C -O -t exon -g \
#     gene_name -a /rumi/shams/genomes/hg38/Homo_sapiens.GRCh38.87.gtf \
#     -o $out $f;
# done

# ############################### samtools ############################################# #
# # Merge .bam files:
# samtools merge -@ 12 bam/all.merged.bam bam/*.bam
# # sort the merged file 
# samtools sort -@12 bam/all.merged.bam > bam/all.merged.srt.bam
# # bam to bet
# bedtools bamtobed -i bam/all.merged.srt.bam > all.merged.srt.bed

# ############################### Systematic ########################################### #
# ref: https://github.com/smithlabcode/piranha
# Piranha -s -b 500 all.merged.srt.bed -o piranha/all.merged.piranha.txt
# # make new ref file 
# awk '{printf("%s\t%d\t%d\t%d\t1\t%s\n",$1,$2,$3, NR,$4)}' piranha/all.merged.piranha.txt > piranha/final_piranha.bed
# awk '{printf("%s\t%s\t%d\t%d\t%s\n",$4,$1,$2+1,$3,$6)}' piranha/final_piranha.bed > piranha/final_piranha.saf
# make new reletive counts
# for f in bam/*_*.bam; do 
#     out=${f/.bam/.piranha.fc}
#     out=${out/bam/piranha}
#     featureCounts -T 12 -O -F SAF -a piranha/final_piranha.saf -o $out $f
# done 

# ############################### Systematic ########################################### #
# for f in bam/*_*.bam; do 
#     out=${f/.bam/.bg}
#     out=${out/bam/bedgraph}
#     # -bga Reporting genome coverage for all positions in BEDGRAPH format.¶
#     bedtools genomecov -bga -ibam $f > $out
# done

# ###############################  meRIPseq  ########################################### #
# https://github.com/goodarzilab/Abe/tree/master/Workflows/meRIPseq/main.R
# ###############################  RHHMM  ############################################## #
# # source: https://github.com/lzcyzm/RHHMM
