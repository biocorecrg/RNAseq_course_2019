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

* Gene Universe
* List of genes selected from the universe

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

### GOSeq

https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf

### GOstats

https://www.bioconductor.org/packages/release/bioc/vignettes/GOstats/inst/doc/GOstatsHyperG.pdf

### Enrichment based on ranked lists of genes

#### GSEA (Gene Set Enrichment Analysis)

**Enrichment score**

Reflects the degree to which a set is overrepresented at the extremes (top or bottom) of the entire ranked list.


