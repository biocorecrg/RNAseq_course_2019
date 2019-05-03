---
layout: page
title: Module 4
navigation: 6
---

# Module 4

## Visualization in the genome browser

### Differences in chromosome nomenclature!

### Conversion to bigwig



## Functional analysis

### Resources

#### Gene Ontology

http://geneontology.org/

#### KEGG pathways

https://www.genome.jp/kegg/

#### Molecular Signatures Database (MSigDB)

http://software.broadinstitute.org/gsea/msigdb/index.jsp


### Enrichment analysis based on gene selection

* Gene Universe: in our example: all genes present in our annotation, on chromosome 10
* List of genes selected from the universe: our selection of genes, give the criteria we defined: padj < 0.01, |log2FoldChange| >= 1

#### Hypergeometric test

Null hypothesis: the list of DEGs is randomly found in the GO.
<br>
Alternative hypothesis: the list of DEGs is over- or under- represented in the GO.
<br>
Universe: all genes in experiment.<br>
Successes: DEGs genes in experiment.<br>
All GO: all genes in GO term.<br>
Successes in GO: DEGs genes in GO.<br>



#### enrichR

http://amp.pharm.mssm.edu/Enrichr/


#### GOSeq

https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf

#### GOstats

https://www.bioconductor.org/packages/release/bioc/vignettes/GOstats/inst/doc/GOstatsHyperG.pdf

#### topGO

https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

library("topGO")
library("org.Hs.eg.db")

     ## GO to Symbol mappings (only the BP ontology is used)
     xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "entrez") # or symbol, ensembl, ...

     allGenes <- unique(unlist(xx)) # retrieve all genes

     myInterestedGenes <- sample(allGenes, 500)

     geneList <- factor(as.integer(allGenes 

     names(geneList) <- allGenes
     
     GOdata <- new("topGOdata",
                   ontology = "BP",
                   allGenes = geneList,
                   nodeSize = 5,
                   annot = annFUN.org, 
                   mapping = "org.Hs.eg.db",
                   ID = "entrez") 


### Enrichment based on ranked lists of genes

#### GSEA (Gene Set Enrichment Analysis)


<img src="images/gsea_presentation.png" width="800" align="middle" />

<img src="images/gsea_msig_banner.png" width="1000" align="middle" />

<img src="images/gsea_msig_sets.png" width="600" align="middle" />




**Enrichment score**

Reflects the degree to which a set is overrepresented at the extremes (top or bottom) of the entire ranked list.


