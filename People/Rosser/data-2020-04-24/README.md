This is how I used [ENCODE Transcription Factor and Histone ChIP-Seq processing pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline2) for current ChIP-Seq dataset. 
## [Installation](https://github.com/ENCODE-DCC/chip-seq-pipeline2#installation)
1. I've cloned the github repo 
```bash
$ cd /rumi/shams/abe/Workflows/
$ git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2
```

2. I made the conda environment as described in this [link](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/install_conda.md). However, I applied uninstallation command before actual installation command to make it work. 
```bash 
$ bash scripts/uninstall_conda_env.sh
$ bash scripts/install_conda_env.sh


3. Activate env
```bash
$ conda activate encode-chip-seq-pipeline
```
3. I downloaded genome database as described [here](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/build_genome_database.md) to build `hg38` genome using *genome builder* command (see [here](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/build_genome_database.md)). 

```bash
$ bash /rumi/shams/abe/Workflows/chip-seq-pipeline2/scripts/download_genome_data.sh hg38 /rumi/shams/abe/genomes/hg38/chip-seq-pipeline2
$ bash /rumi/shams/abe/Workflows/chip-seq-pipeline2/scripts/build_genome_data.sh hg38 /rumi/shams/abe/genomes/hg38/chip-seq-pipeline2
```

4. After running builder command, there is a TSV file which is required to define input JSON file. 
```
    "chip.genome_tsv" : "/rumi/shams/abe/genomes/hg38/chip-seq-pipeline2/hg38.tsv"
```

* The genome files and indecies are `/rumi/shams/abe/genomes/hg38/chip-seq-pipeline2/` which is ready for other people to use as well. 

## Running the pipeline
1. Make sure to activate the conda environment:
```bash
$ conda activate encode-chip-seq-pipeline 
```

2. Prepare inputs:
All raw fastq files stored at a child `fastq` directory at `/rumi/shams/abe/People/Rosser/data-2020-04-24/`. Then, I made four json files in four seprate directories (RT112-H3K27ac, RT112-H3K4me3, UMUC3-H3K27ac, UMUC3-H3K4me3) as described [here](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/input_short.md) using the [template](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/example_input_json/template.json). 

3. I just follow steps 3, 5 and 6 to [install caper](https://github.com/ENCODE-DCC/caper#installation). `rumi` server is a `local` platform. 

```bash
$ caper init local
```

See this [link](https://github.com/ENCODE-DCC/caper#running-pipelines-on-general-computers) for more information about using capper.

4. I'm switching to each four dirctories to run the pipeline seperatly. For example: 
```bash
$ cd /rumi/shams/abe/People/Rosser/data-2020-04-24/RT112-H3K27ac
$ nohup caper run ~/Workflows/chip-seq-pipeline2/chip.wdl -i input.jason > caper.log
```

## [Organize outputs](https://github.com/ENCODE-DCC/chip-seq-pipeline2#how-to-organize-outputs)
1. After the complete run of the pipeline, there is a `JSON` file at `chip/some-hash-code/meta.json` path. `croo` command organize outputs using this file. Again, I do this step four times for each study. For example:
```bash
$ cd /rumi/shams/abe/People/Rosser/data-2020-04-24/RT112-H3K27ac
$ croo chip/some-hash-code/metadata.json
```

2. Merge qc reports from different studies:
At the same conda envirnment `pip install qc2tsv` install (qc2tsv)[https://github.com/ENCODE-DCC/qc2tsv]. Then, it's ready to use for merging qc reports from different studies in a unique spreadsheet file.
```bash 
$ cd /rumi/shams/abe/People/Rosser/data-2020-04-24/
$ qc2tsv RT112-H3K27ac/qc/qc.json RT112-H3K4me3/qc/qc.json UMUC3-H3K27ac/qc/qc.json UMUC3-H3K4me3/qc/qc.json > spreadsheet.tsv
```

~Abe, May 2020
