# ExpressionAnalysisPipeline
Pipeline that receives input of multiple fastq files and returns a tsv file with expression level, Entrez ID, GO term etc for each gene.

## Require Tools
- Python3
- Trimmomatic
- HISAT2
- samtools
- stringtie
- mygene(python package)

## Pipeline Flow Diagram
![flow.png](https://github.com/suecharo/ExpressionAnalysisPipeline/raw/master/flow.png)

## Input Files
### Fastq Files
- Multiple single-end · paired-end files can be received at once.
- Avaliable extensions are `.fq` , `.fastq` , `.fq.gz` .
- In the case of paired-end, connect two file paths with a comma.

### Annotation Files
- You have to specify either Fasta File or build completed HISAT2 index files.
- You must enter the gtf / gff file.

### Trimmomatic Adapter Path.
- You must set conf.yml in the same directory and describe the paths of Trimmomatic path and adapter path in it.

``` conf.yml
TRIMMOMATIC_PATH : /usr/local/bioinfo_tools/Trimmomatic-0.36/trimmomatic-0.36.jar
TRIM_ADAPTER_SE : /usr/local/bioinfo_tools/Trimmomatic-0.36/adapters/single_end.fa
TRIM_ADAPTER_PE : /usr/local/bioinfo_tools/Trimmomatic-0.36/adapters/paired_end.fa
```

## Usage
### Help Message

```
$ python3 expression_analysis.py --help
usage: expression_analysis [-h] [-a FASTA | -i INDEX] -g GTF [-o OUTPUT]
                           [-p CPU]
                           fastq [fastq ...]

Expression analysis pipeline created by suecharo.

positional arguments:
  fastq       Input all fastq file.(.fastq or .fa or .fa.gz) (single-
              end=SRR000001.fq, paired-end=SRR000002_1.fq,SRR000002_2.fq)

optional arguments:
  -h, --help  show this help message and exit
  -a FASTA    Input fasta file of .fa or .fasta.
  -i INDEX    Input HISAT2 index path.(e.g. ./hisat2-index/hg19)
  -g GTF      Input gtf file of hg19 or mm9.
  -o OUTPUT   Enter output dir.(default=./output)
  -p CPU      Input cpu num.(default=4)
```

### Execution Command Example

``` bash
$ python3 expression_analysis.py -i /home/suecharo/data/expression_analysis/hisat2_index/hg19 -g /home/suecharo/data/expression_analysis/hg19.gtf -o /home/suecharo/analysis/expression_analysis/test -p 20 /home/suecharo/data/expression_analysis/data/SRR951071/SRR951071.fastq /home/suecharo/data/expression_analysis/data/SRR951072/SRR951072_1.fastq,/home/suecharo/data/expression_analysis/data/SRR951072/SRR951072_2.fastq
```

### Differences Between Each Script
- expression_analysis.py
    - Process all isoforms.
- expression_analysis_longest.py
    - Process the longest isoform.

### Column of Output Data
- RefSeq_ID
- each_FPKM
- each_TPM
- each_Coverage
- Gene_ID
- Chromosome
- Start
- End
- Width
- Strand
- Entrez_Gene_ID
- Description
- Cellular_Component_Accession
- Cellular_Component_Name
- Molecular_Function_Accession
- Molecular_Function_Name
- Biological_Process_Accession
- Biological_Process_Name

### Output Example
![output_1.png](https://github.com/suecharo/ExpressionAnalysisPipeline/raw/master/output_1.png)

![otuput_2.png](https://github.com/suecharo/ExpressionAnalysisPipeline/raw/master/output_2.png)
