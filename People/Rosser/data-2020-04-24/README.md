[**ENCODE Transcription Factor and Histone ChIP-Seq processing pipeline**](https://github.com/ENCODE-DCC/chip-seq-pipeline2)

This is how I used above pipeline for current ChIP-Seq dataset. 
1. [Pipeline installation](https://github.com/ENCODE-DCC/chip-seq-pipeline2#installation)
I've cloned the github repo. Then, I made the conda enviernment as mentioned in this [link](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/install_conda.md). I applied uninstallation command before actual installation command. 
Finally, I downloaded genome database as described [here](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/build_genome_database.md) to build `hg38` genome. It's in the rumi server at: `/rumi/shams/abe/genomes/hg38/chip-seq-pipeline2/`

2. Prepare inputs 
All raw fastq files stored at a `fastq` directory at `/rumi/shams/abe/People/Rosser/data-2020-04-24/`. 
Then, I made four json files in four seprate directories (RT112-H3K27ac, RT112-H3K4me3, UMUC3-H3K27ac, UMUC3-H3K4me3) as described [here](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/input_short.md) using the [template](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/example_input_json/template.json). 
I'm following this [check list](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/input_short.md#checklist). 

3. Run pipeline
First, I activate the envirnment: `conda activate encode-chip-seq-pipeline`. 
Then, I just follow steps 3, 5 and 6 to [install caper](https://github.com/ENCODE-DCC/caper#installation). (rumi is a `local` platform). 
See this [link](https://github.com/ENCODE-DCC/caper#running-pipelines-on-general-computers) for using capper. 

I'm following this [check list](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/input_short.md#checklist)

4. [How to organize outputs](https://github.com/ENCODE-DCC/chip-seq-pipeline2#how-to-organize-outputs)
After the complete run of the pipeline, there is a `JSON` file at `study-name/chip/some-hash-code/meta.json`. You use `croo` command with this file to organize outputs. 

`cd study-name`
`croo chip/some-hash-code/metadata.json`

5. Merge qc reports from different studies 
At the same conda envirnment `pip install qc2tsv` to merge all qc reports in a unique spreadsheet file. So, I used this command at the parent directory to the fastq files and all outputs:
`qc2tsv RT112-H3K27ac/qc/qc.json RT112-H3K4me3/qc/qc.json UMUC3-H3K27ac/qc/qc.json UMUC3-H3K4me3/qc/qc.json > spreadsheet.tsv`

~Abe, May 2020
