---
layout: page
title: Module 3
navigation: 5
---

# Module 3

## Experimental design

### Number of replicates

### Depth of sequencing



## Differential expression analysis

Differential expression analysis means to perform statistical analysis to try and discover changes in expression levels (genes, transcripts) between experimental groups.<br>



### DESeq2

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is an R/Bioconductor implemented method to detect differentially expression genes.<br>
Bioconductor is an R project and repository that provides a set of packages and methods for omics data analysis.<br>
It is based on the negative binomial distribution.<br>
Describes the distribution of draws with replacement. The number of failures is fixed.
<br>
The package DESeq2 provides methods to test for differential expression by use of negative binomial generalized linear models.

This DESeq2 tutorial is widely inspired from the [RNA-seq workflow](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html) developped by the authors of the tool.


* Raw count matrices

DESeq2 expects a matrix of integer values. The value at the i-th row and j-th column tells how many reads have been assigned to gene i in sample j.

| gene name | sample A | sample B | sample C | sample D | sample E | sample F |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| gene 1 | 0 | 10 | 5 | 4 | 1 | 2 |
| gene 2 | 123 | 189 | 76 | 38 | 132 | 139 |
| gene 3 | 20 | 25 | 35 | 43 | 27 | 32 |
| gene 4 | 67 | 72 | 32 | 79 | 62 | 59 |
| gene 5 | 230 | 241 | 229 | 120 | 99 | 103 |

* Sample sheet

The sample sheet describes the samples characteristics: treatment, knock-out / wild type, replicates etc.

| SampleName | CellType | Group |
| :---: | :---: | :---: |
| sample A | HeLa | WT |
| sample B | HeLa | WT |
| sample C | HeLa | WT |
| sample D | HeLa | KO |
| sample E | HeLa | KO |
| sample F | HeLa | KO |

* Load data from STAR

```
mkdir counts_star
```

The "ReadsPerGene" output files of STAR (from option --quantMode GeneCounts) contain 4 columns that correspond to different counts / read overlap according to the protocol's strandedness:
* column 1: gene ID
* column 2: counts for unstranded RNA-seq.
* column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
* column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse): the most common protocol.

for i in `ls -d */`
do echo $i
# retrieve the gene ID and 2nd read strand count information corresponding columns
cut -f1,4 $i/*ReadsPerGene.out.tab | grep -v "_" > counts_star/`basename $i/*PerGene.out.tab ReadsPerGene.out.tab`_counts.txt
done

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ group)

se <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                        directory = countsdir,
                        design = ~ group)


* Load counts data from SALMON

files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- samples$run
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

* 

se1 <- DESeq(se)

# filter out rows which have sum 0 or 1 across all samples
se2 <- se1[rowSums(counts(se1)) > 1, ] # 25393 versus 50600 the initial 


* Samples correlation

Sample-to-sample distances

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

* Principal Component Analysis

Reduction of dimensionality to be able to retrieve main differences between samples

plotPCA(vsd, intgroup = "group")


* Running differential expression analysis

comp_tmp <- results(se3, contrast=c("group", "KO", "WT"), cooksCutoff=0.20)
# KO vs WT

The column log2FoldChange is the effect size estimate. It tells us how much the gene’s expression seems to have changed due to treatment with dexamethasone in comparison to untreated samples. This value is reported on a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene’s expression is increased by a multiplicative factor of 21.5≈2.82.

DESeq2 performs for each gene a hypothesis test to see whether evidence is sufficient to decide against the null hypothesis that there is zero effect of the treatment on the gene and that the observed difference between treatment and control was merely caused by experimental variability (i.e., the type of variability that you can expect between different samples in the same treatment group). As usual in statistics, the result of this test is reported as a p value, and it is found in the column pvalue. Remember that a p value indicates the probability that a fold change as strong as the observed one, or even stronger, would be seen under the situation described by the null hypothesis.

* Running the R script from the command line

Rscript deseq2_star.R sampletable_star.txt
Rscript deseq2_salmon.R sampletable_salmon.txt



### Online tool

This [online tool](http://52.90.192.24:3838/rnaseq2g/) provides a way to process differential expression analysis using some of the popular tools in the field (among whichDESeq2, edgeR, limma), starting from raw counts, using an UI (User Interface).

* Prepare the matrix of raw counts

The first column contains, in our case, the gene names.
The remaining columns contain the expression of each gene in each sample (one column per sample).

* Prepare the sample sheet

Sample names must match column names in matrix.<br>
Add one column that corresponds to the experimental groups the samples belong to


* Choose the control group name and the case group name

In our case:<br>
Control: WT<br>
Case: KO<br>

* Choose which samples belong to which experimental group

Control samples: <br>
Case samples: <br>

* Paired test

* Choose DE method(s)

Let's try DESeq2 and edgeR

* Submit DE analysis

* Go to "Results" tab



## Gene selection

* padj (p-value corrected for multiple testing)
* log2FC (log2 Fold Change)


## Visualization of differential expression

### Volcano plots

```
# code for volcano plot
```

Run the R script:<br>
Rscript volcano_plot.R


### Individual genes






