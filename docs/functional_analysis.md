---
layout: page
title: Functional analysis
navigation: 29
---

# Functional analysis

## Data bases

### Gene Ontology

<img src="images/GO_logo.png" width="200" align="middle" />

The [Gene Ontology (GO)](http://geneontology.org/) describes our knowledge of the biological domain with respect to three aspects:

| GO domains / root terms | Description |
| :----: | :----: |
| Molecular Function | Molecular-level activities performed by gene products. e.g. **catalysis**, **binding**. |
| Biological Process | Larger processes accomplished by multiple molecular activities. e.g. **apoptosis**, **DNA repair**. |
| Cellular Component | The locations where a gene product performs a function. e.g. **cell membrane**, **ribosome**. |

Example of GO annotation: the gene product "cytochrome c" can be described by the **molecular function** *oxidoreductase activity*, the **biological process** *oxidative phosphorylation*, and the **cellular component** *mitochondrial matrix*.
<br><br>
The structure of GO can be described as a graph: each GO term is a **node**, each **edge** represents the relationships between the nodes. For example:

<img src="images/GO_example_graph.png" width="500" align="middle" />

GO:0019319 (hexose biosynthetic process) is part of GO:0019318 (hexose metabolic process) and also part of  GO:0046364 (monosaccharide biosynthetic process). They all share common **parent nodes**, for example GO:0008152 (metabolic process), and eventually a **root node** that is here **biological process**.

### KEGG pathways

The [Kyoto Encyclopedia of Genes and Genomes (KEGG)](https://www.genome.jp/kegg/) is a database for understanding high-level functions and utilities of the biological system. <br>

It provides comprehensible manually-drawn pathways representing biological processes or disease-specific pathways.<br>
Example of the [*Homo sapiens* melanoma pathway](https://www.genome.jp/dbget-bin/www_bget?hsa05218):

<img src="images/kegg_hsa05218.png" width="700" align="middle" />


### Molecular Signatures Database (MSigDB)

<img src="images/gsea_msig_banner.png" width="1000" align="middle" />

The [Molecular Signatures Database (MSigDB)](http://software.broadinstitute.org/gsea/msigdb/index.jsp) is a collection of 17810 annotated gene sets (as of May 2019) created to be used with the GSEA software (but not only). <br>

It is divided into [8 major collections](http://software.broadinstitute.org/gsea/msigdb/collections.jsp) (that include the previously explained GOs and KEGG pathways):

<img src="images/gsea_msig_sets.png" width="300" align="middle" />


## Enrichment analysis based on gene selection

Tools based on a user-selection of genes usually require 2 inputs:

* Gene Universe: in our example: all genes used in our analysis (after filtering out low counts in our case).
* List of genes selected from the universe: our selection of genes, give the criteria we previously used: **padj < 0.05**, **&#124;log2FoldChange&#124; >= 0.5**.

### GO / Panther tool

The [main page of GO](http://geneontology.org/) provides a tool to test the enrichment of gene ontologies or Panther/Reactome pathways in pre-selected gene lists.
<br>
The tool needs a selection of differentially expressed genes (supported IDs are: gene symbols, ENSEMBL IDs, HUGO IDs, UniGene, ..) and a gene universe.
<br>

Prepare files using the ENSEMBL IDs:

```{bash}
# Extract all gene IDs used in our analysis and convert from Gencode (e.g ENSG00000057657.16) to ENSEMBL (e.g. ENSG00000057657) IDs
cut -f1 deseq2_results.txt | sed '1d' | sed 's/\..//g' > deseq2_universe_ensemblIDs.txt

# Convert from Gencode to ENSEMBL IDs from selected gene list
sed 's/\..//g' deseq2_results_padj0.05_log2fc0.5_IDs.txt > deseq2_results_padj0.05_log2fc0.5_ensemblIDs.txt
```
Paste our selection, and select **biological process** and **Homo sapiens**(file deseq2_results_padj0.05_log2fc0.5_IDs.txt):

<img src="images/GO_tool_interface.png" width="500" align="middle" />

**Launch** !

<img src="images/GO_tool_input1.png" width="800" align="middle" />

**Analyzed List** is what we just uploaded (*deseq2_results_padj0.05_log2fc0.5_ensemblIDs.txt*).
<br>
In **Reference List**, we need to upload a file containg the **universe** (*deseq2_universe_ensemblIDs.txt*):  *Change -> Browse -> (select deseq2_universe_ensemblIDs.txt) -> Upload list*
<br>

* Launch analysis
<img src="images/GO_tool_results_ensembl.png" width="800" align="middle" />
* Try the same analysis using the **gene symbols** instead of ENSEMBL IDs
```{bash}
# Get universe with gene symbols (we already have the gene selection in deseq2_results_padj0.05_log2fc0.5_symbols.txt)
cut -f2 deseq2_results.txt | sed '1d' > deseq2_universe_symbols.txt
```
* **Launch** !
<img src="images/GO_tool_results_symbols.png" width="800" align="middle" />



### enrichR

[EnrichR](http://amp.pharm.mssm.edu/Enrichr/) is a gene-list enrichment tool developped at the Icahn Schoold of Medicine (Mount Sinai).

<img src="images/enrichr_interface.png" width="500" align="middle" />

It does not require the input of a gene universe: only a selection of genes or a BED file.

<img src="images/enrichr_paper1.jpg" width="500" align="middle" />

The default EnrichR interface works for *Homo sapiens* and *Mus musculus*.<br>
However, EnrichR also provides a [set of tools](https://amp.pharm.mssm.edu/modEnrichr/) for ortholog conversion and enrichment analysis of more organisms:

<img src="images/enrichr_interface2.png" width="600" align="middle" />

In the [main page](http://amp.pharm.mssm.edu/Enrichr/), paste our list of selected **gene symbols** (*deseq2_results_padj0.05_log2fc0.5_symbols.txt*) and **Submit** !

<img src="images/enrichr_results_all.png" width="600" align="middle" />

KEGG Human pathway **bar graph** vizualization:
<img src="images/enrichr_results_bar.png" width="600" align="middle" />

KEGG Human pathway **table** vizualization:
<img src="images/enrichr_results_table.png" width="600" align="middle" />

KEGG Human pathway **clustergram** vizualization:
<img src="images/enrichr_results_clustergram.png" width="600" align="middle" />

For **Cell Types**, you can also visualize networks, for example **Human gene Atlas**:
<img src="images/enrichr_results_network.png" width="600" align="middle" />

You can also export some graphs as PNG, JPEG or SVG.


## Enrichment based on ranked lists of genes using GSEA

### GSEA (Gene Set Enrichment Analysis)

<img src="images/gsea_presentation.png" width="600" align="middle" />

[GSEA](http://software.broadinstitute.org/gsea/) is available as a Java-based tool.

#### Algorithm

GSEA doesn't require a threshold: the whole set of genes is considered.

<img src="images/gsea_paper.jpg" width="800" align="middle" />

GSEA checks whether a particular gene set (for example, a gene ontology) is **randomly distributed** across a list of **ranked genes**.
<br>
The algorithm consists of 3 key elements:

1. **Calculation of the Enrichment Score**
The Enrichment Score (ES) reflects the degree to which a gene set is overrepresented at the extremes (top or bottom) of the entire ranked gene list.
2. **Estimation of Significance Level of ES** 
The statistical significant (nominal p-value) of the **Enrichment Score (ES)** is estimated by using an empirical phenotype-based permutation test procedure. The **Normalized Enrichment Score (NES)** is obtained by normalizing the ES for each gene set to account for the size of the set.
3. **Adjustment for Multiple Hypothesis Testing**
Calculation of the FDR ti control the proportion of falses positives.

<img src="images/gsea_explained.gif" width="900" align="middle" />

See the [GSEA Paper](https://www.ncbi.nlm.nih.gov/pubmed/16199517) for more details on the algorithm.

The main GSEA algorithm requires 3 inputs:
* Gene expression data
* Phenotype labels
* Gene sets

#### Gene expression data in TXT format

The input should be normalized read counts filtered out for low counts (-> we created it in the DESeq2 tutorial -> *normalized_counts.txt* !).
<br>
The first column contains the gene ID (HUGO symbols for *Homo sapiens*).<br>
The second column contains any description or symbol, and will be ignoreed by the algorithm.<br>
The remaining columns contains normalized expressions: one column per sample.

| NAME | DESCRIPTION | A549_0_1 | A549_0_2 | A549_0_3 | A549_25_1 | A549_25_2 | A549_25_3 |
| DKK1 | NA| 0 | 0 | 0 | 0 | 0 | 0 |
| HGT | NA | 0 | 0 | 0 | 0 | 0 | 0 |

<b>Exercise</b>
<br>
Adjust the file **normalized_counts.txt** so the first column is the gene symbol, the second is the gene ID (or anything else), and the remaining ones are the expression columns. You can save that new file as **gsea_normalized_counts.txt**.
<br>

```{bash}
awk 'BEGIN{OFS="\t"}{print $2,$1,$3,$4,$5,$6,$7,$8}' normalized_counts.txt > gsea_normalized_counts.txt

```


#### Phenotype labels in CLS format

A phenotype label file defines phenotype labels (experimental groups) and assigns those labels to the samples in the corresponding expression data file.

<img src="images/gsea_phenotypes.png" width="700" align="middle" />

Let's create it for our experiment:

| 6 | 2 | 1 |  |  |  |
| # | t0 | t25 |  |  |  |
| t0 | t0 | t0 | t25 | t25 | t25 |

**NOTE**: the first label used is assigned to the first class named on the second line; the second unique label is assigned to the second class named; and so on. 
<br> 
So the phenotype file could also be:

| 6 | 2 | 1 |  |  |  |
| # | t0 | t25 |  |  |  |
| 0 | 0 | 0 | 1 | 1 | 1 |

The first label **t0** in the second line is associated to the first label **0** on the third line.

<b>Exercise</b>
<br>
Create the phenotype labels file and save it as **gsea_phenotypes.cls**.


Reflects the degree to which a set is overrepresented at the extremes (top or bottom) of the entire ranked list.


#### Download and run GSEA

##### Download Java application:

Enter the [download page](http://software.broadinstitute.org/gsea/login.jsp), enter your Email and **login**:
<img src="images/gsea_downloads.png" width="800" align="middle" />
<br>
Click on **download gsea-3.0.jar** link and save file locally to your home directory.


##### Launch the GSEA application

GSEA is Java-based. Launch it from a terminal window:

```{bash}
$RUN java -Xmx1024m -jar gsea-3.0.jar
```

<img src="images/gsea_interface.png" width="800" align="middle" />

In the *Steps in GSEA analysis* (upper left corner):

* Go to **Load data**: select **gsea_normalized_counts.txt** and **gsea_phenotypes.cls** and load.

<img src="images/gsea_load_files.png" width="800" align="middle" />

* Go to **Run GSEA**

<img src="images/gsea_parameters_commented.png" width="800" align="middle" />

* Results: **index.html**

<img src="images/gsea_results_index_commented.png" width="800" align="middle" />

* **enrichments results in html**

<img src="images/gsea_results_stats.png" width="800" align="middle" />

* Details

<img src="images/gsea_results_details1.png" width="800" align="middle" />

<img src="images/gsea_results_details2.png" width="600" align="middle" />

<img src="images/gsea_results_details3.png" width="500" align="middle" />

------------

PENDING FINISHING
