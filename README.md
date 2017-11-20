#RNA-seq analysis 
11.19.17
Jesse Dabney

###Required programs

* FastQC
* samtools
* bioawk
* Trim Galore!
   * requires cutadapt and FastQC
* HISAT2
* HTSeq
* R  
   * tidyverse
   * edgeR


###download data

I'm using toy data from a study involving commercially available human RNA from cancer cell lines and brain tissue. Each set contains 3 technical replicates and has been pre-filtered for reads mapping to chromosome 22.

Sequence data, fastq

```bash
mkdir -p ./trial1/data && cd $_
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar 
```

reference genome, fasta

```bash
wget http://genomedata.org/fasta/GRCh38/chr22_with_ERCC92.fa
samtools faidx chr22_with_ERCC92.fa
```

gene information, gtf

```bash
wget http://genomedata.org/annotations/GRCh38/chr22_with_ERCC92.gtf
```

check number of reads in each file

```bash
for i in *.fastq.gz; do echo $i; bioawk -cfastx 'END{print NR}' $i; done
```


###pre-trimming and alignment QC

Here I'm running FastQC on the downloaded fastq files. FastQC gives some useful stats that can be used to assess the quality of the reads, including base and sequence quality scores, GC content and overrepresented sequences (e.g. adapters or contamination).

```bash
for i in *.fastq.gz; do fastqc -o ./fastqc ../data/$i; done
```

This produces html files, which I visually inspected. It also outputs files that can be parsed as part of a pipeline, so that samples failing this step can be removed.


###Adaptor Trimming

I'm using Trim Galore! here, based on prior experience with this program. It will trim off low quality bases and adapter sequences, and also run FastQC again on the trimmed reads.

A brief explanation of the options:

  * -q 15: removes low quality bases from the ends of reads
  * --fastqc_args -t 4: run fastqc with 4 cores
  * --illumina: look for illumina adapter sequences. The libraries were produced with TruSeq RNA kits.
  * --retain_unpaired: generate separate file with unpaired reads
  * --paired: specify that the data is paired end

```bash
trim_galore -q 15 --fastqc_args '-t 4' --illumina --retain_unpaired --paired ../data/HBR_Rep3*read1*.gz ../data/HBR_Rep3*read2*.gz
trim_galore -q 15 --fastqc_args '-t 4' --illumina --retain_unpaired --paired ../data/UHR_Rep3*read1*.gz ../data/UHR_Rep3*read2*.gz
```

###Alignment

I'm using HISAT2 for alignment. Some other options were Tophat and STAR, but it seems that HISAT2 is an improvement to Tophat, and STAR was too memory intensive to run on my laptop.

HISAT2 can work in a splicesite aware manner, so first I need to get the splicesites, using the gtf file. This script outputs the start and end position of the site, and which strand it's on.

```bash
hisat2_extract_splice_sites.py ./chr22_with_ERCC92.gtf > chr22_splicesites.tsv
```

And then I need to get a list of exons. This will have the same format as the splicesite file.

```bash
hisat2_extract_exons.py chr22_with_ERCC92.gtf > chr22_exons.tsv
```

HISAT2 also needs an index.

```bash
hisat2-build -p 4 --ss ./chr22_splicesites.tsv --exon ./chr22_exons.tsv chr22_with_ERCC92.fa ./chr22_with_ERCC92
```

Now I run the alignment. In the interest of time I ran these individually, but this can easily be looped. Output is a sam file. An example command, with options:

* -p 4: run on 4 cores
* --rg*: specify strings for fields in bamfile
* -x: directory with index
* --dta: transcript alignment
* --rna-strandedness: RF since libraries made with strand specific kit

```bash
hisat2 -p 4 --rg-id=HBR_Rep3 --rg SM:HBR --rg LB:HBR_Rep3_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x ../refs/chr22_with_ERCC92 --dta --rna-strandness RF -1 ../trimming/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1_val_1.fq.gz -2 ../trimming/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2_val_2.fq.gz -S ./HBR_Rep3.sam
```

Sort sam file and convert to bam

```bash
for i in *.sam; do samtools sort -@ 4 -o ${i/sam/bam} $i; done
```

Merge bamfiles

```bash
samtools merge HBR.bam HBR_Rep1.bam HBR_Rep2.bam HBR_Rep3.bam
```

Index bamfile

```bash
samtools index HBR.bam
```

###Inspect alignment

It's good to have a look at the alignments, for instance looking at the samflags for read fates, or generating summary stats. Some things that can be done here:

* make mpileup files
* samtools view
* samtools tview
* samtools flagstat

###Get expression profile with HTSeq (counting reads mapped to genes)

To prepare for differential gene expression analysis, I am generating counts of reads that map to exons. 

```bash
for i in HBR_R*.bam; do htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $i ../refs/chr22_with_ERCC92.gtf > ../htseq_counts/${i/.bam/_gene.tsv}; done
```

And then merge results to one file

```bash
join UHR_Rep1_gene.tsv UHR_Rep2_gene.tsv | join - UHR_Rep3_gene.tsv | join - HBR_Rep1_gene.tsv | join - HBR_Rep2_gene.tsv | join - HBR_Rep3_gene.tsv > gene_read_counts_table_all.tsv
```


###DGE analysis with edgeR

I chose edgeR for 2 reasons. One, I'm fairly comfortable with R so I figured it would show some proficiency in something other than bash, and two, it was mentioned in the job description.

Get started with some upfront R stuff:

Install the necessary libraries and set working directory:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
library(tidyverse)
library(magrittr)

workdir <- "~/Desktop/work_projects/rna_seq_trial/trial1/htseq_counts/"
setwd(workdir)
```

Then read in the data generated above:

```r
gene.names <- read_tsv(file = "./ENSG_ID2Name.txt", col_names = F)
data <- read_tsv(file = "./gene_read_counts_table_all_final.tsv")
```

Generate some isolated lists for the analysis (sample, geneIDs, gene names).

```r
labs <- factor( c( rep("UHR",3), rep("HBR",3) ))
genes <- data %>% select(GeneID) %>% pull()
g.names <- gene.names %>% filter(X1 %in% genes) %>% pull(X2)
```

The data needs to be converted from tibble to matrix for edgeR.

```r
data.m <- data %>% select(2:7) %>% as.matrix
rownames(data.m) <- data %>% pull(GeneID)
```

And now calculate DGE with edgeR. 

```r
y <- DGEList(counts = data.m, genes = g.names, group = labs)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
topTags(et)
de <- decideTestsDGE(et, p=.05)
detags <- rownames(y)[as.logical(de)]
```

And then combine the names of genes with significantly different levels of expression (in either direction) with p-values and fold-chages. This file is called DE_genes.txt and is located in this repo.

```r
mat <- cbind(
  genes,g.names,
  sprintf('%0.3f',log10(et$table$PValue)),
  sprintf('%0.3f',et$table$logFC)
)[as.logical(de),]

colnames(mat) <- c("Gene", "Gene_Name", "Log10_Pvalue", "Log_fold_change")

o <- order(et$table$logFC[as.logical(de)],decreasing=TRUE)
mat <- mat[o,]
write.table(mat, file="DE_genes.txt", quote=FALSE, row.names=FALSE, sep="\t")
```

###GO Enrichment analysis

With the list of genes, I want to know what the functional difference is between the two sample sets. This can be accomplished with GO enrichment analysis.

It looked like edgeR has functionality to do this, but I was short on time so opted to use an online plug-and-go option: http://www.geneontology.org/page/go-enrichment-analysis. I generated two files, one looking at biological process categories, and another looking at cellular component categories. The results indicate enrichment for genes involved in protein transport and localization. 

Since this analysis is based only on chromosome 22 genes, the results highly informative. A naive biological interpretation would be that different demands on cancer and brain cells lead to differences in gene expression levels on chromosome 22, and that these differences lie in genes involved in protein movement. This may have something to do with signalling pathways in brain tissue, or cell adhesion in cancer.






