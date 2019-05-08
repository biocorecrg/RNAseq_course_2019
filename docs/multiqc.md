---
layout: page
title: MultiQC report
navigation: 15
---

# Combining reports
At this point, we can summarize all the work done by using the tool [**MultiQC**](https://multiqc.info/). 
MultiQC aggregates outputs from many bioinformatics tools across many samples into a single report by searching a given directory for analysis logs and compiling a HTML report. 

First, link our mapping results to the directory QC:
```{bash}
cd QC/
ln -s ../alignments/* .
```

Then run multiqc on the directory QC to combine all reports:

```{bash}
$RUN multiqc .
[INFO   ]         multiqc : This is MultiQC v1.7 (7d02b24)
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching 'QC/'
Searching 70 files..  [####################################]  100% 
...
```

and visualize the final report in the browser:

```{bash}
firefox multiqc_report.html
```



<img src="images/multiqc.png"  align="middle" />

<br/>
