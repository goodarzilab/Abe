# ############################### Trimming ############################################# #
# cd ~/People/Judd/fastq
# for f in *_L002_R1_001.fastq.gz; do  
# out=${f/\_S[1-9]*/.trim.fastq.gz}; 
# cutadapt -j 12 -q 15 -m 20 -a NNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o ../trim/$out $f; 
# done
# ############################### Alignment ############################################ #
# STAR --genomeLoad LoadAndExit --genomeDir /rumi/shams/abe/genomes/hg38/
# for f in trim/*.fastq.gz;do
# out=${f/.trim.fastq.gz/_};
# out=${out/trim/bam};
# STAR --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --runThreadN 16 --genomeDir /rumi/shams/abe/genomes/hg38/ --readFilesIn $f --outFileNamePrefix $out --outReadsUnmapped Fastx;
# done
# STAR --genomeLoad Remove --genomeDir /rumi/shams/abe/genomes/hg38/

# ############################### all files renamed #################################### #

# # ############################### Alignment: HIV genome ############################## #
for f in fastq_unmapped/*; do
	o=${f/.fastq/.bam};
 	o=${o/fastq_unmapped/bam};
 	bowtie2 --sensitive -N 1 -x HIV/HIV -U $f | samtools sort -o $o;
done

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

# ###############################  .bed to .fa ######################################### #
# use bioconda ucsc-twobittofa
# conda activate bedtools
# # remove duplicates from bed file
# grep -v '#' peak.bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,"1",$6}' > peak.c.bed
# twoBitToFa -noMask -bed=peak.drop_dup.bed /rumi/shams/genomes/hg38/hg38.2bit peak.fa

# ###############################  .bam to .bigwig ##################################### #
# use bioconda deeptools
# for b in bam/*.bam; do
#   bw=${b/.bam/.bw};
#   bw=${bw/bam/bigwig};
#   bamCoverage -b $b -o $bw;
# done
# ###############################  prep_seqs_for_teiser_run ############################ #
# perl /flash/hani/bin/Tools/prep_seqs_for_teiser_run.pl peak.fa peaks
# ###############################  Run fire ############################################ #
#perl $FIREDIR/fire.pl --expfile=~/People/Judd/exomepeak/hg38/control/peaks_teiser.txt \
#                      --exptype=discrete --fastafile_rna=~/People/Judd/exomepeak/hg38/control/peaks_teiser.fa \
#                      --nodups=1 --dodna=0 --dodnarna=0 --species=human
#
# perl $FIREDIR/fire.pl --expfile=peaks_teiser.txt \
#                       --exptype=discrete --fastafile_rna=peaks_teiser.fa \
#                       --nodups=1 --dodna=0 --dodnarna=0 \
#                       --species=human --doskipdiscovery=1 --motiffile_rna=<FILE>
