# MGSAT semi-automated differential abundance analysis of omics datasets #

![Within-between dissimilarity plot example](http://andreyto.bitbucket.org/mgsat/example_project_01/plots/10584217153b.png "Within-between dissimilarity plot example")

##Overview##
MGSAT is written in R. It applies several types of statistical tests, normalizations and plotting routines to the abundance count matrices that are typically the output of annotating (meta)omics datasets, and generates a [structured HTML report](http://andreyto.bitbucket.org/mgsat/example_project_01/0-report.html) that, in addition to results, shows method parameters and versions of the external packages. 

The user has fine-grained control over types of tests, parameters, and a description of a study design through a named list data structure that is provided as input to the top-level routine of the package.

MGSAT has being used to analyse datasets from 16S gene environmental sequencing surveys, proteomics mass spectrometry and meta-genomics whole-genome shotgun functional annotation in several studies where the goal was to associate the state of the human microbiome or human proteome with a specific disease condition such as leukemia, diabetes or Respiratory syncytial virus (RSV).

For 16S surveys, MGSAT has routines that load output count matrices generated through [Mothur-based Standard Operating Procedure](http://www.mothur.org/wiki/MiSeq_SOP) or compatible automation wrappers such as [YAP](https://github.com/andreyto/YAP). The analysis can be done at all levels of annotated taxonomy or for the operational taxonomic units (OTU).

MGSAT targets study designs where the goal is to associate annotated omics data with clinical metadata variables. Between-, within-subject and mixed designs can be addressed. The analysis includes alpha and beta diversity, richness estimates, both abundance based and incidence based, tests for difference of diversity and richness between groups, PermANOVA test for abundance profile dissimilarities for either independent or paired samples, comparison of paired vs unpaired profile dissimilarities, stability ranking analysis based on rank sum and signed rank Wilcoxon tests, stability selection analysis based on elastic net classifier, [DESeq2](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) differential abundance analysis based on a negative binomial model of abundance counts, and a binomial mixed model analysis where overdispersion is accounted for by a subject-specific random effect. Optional data normalizations includes simple proportions, inverse hyperbolic sine, Box-Cox, Aitchison and DESeq2 transforms.

Plots include abundance profile plots in different representations split by grouping variables or shown as trends along continuous metadata variables, as well as clustered heatmaps with overlayed panels showing metadata variables.

##Sample Analysis Code and Output##
The easiest way to start using this package is to study the inputs and outputs of an example analysis project.

The project is a made-up 16S sequencing study that looked at the association of a gut microbiome with some diet regimen. Samples were taken at several longitudinal points (called "visits" here), each visit coming after the next phase of diet regimen. Visit 1 was before the start of treatment. At Visit 1 (and a few at Visit 2), samples were also taken from healthy control individuals, each of which was matched with a corresponding patient. Some patients received a separate incidental drug treatment before coming for visit 1, which was reflected in a corresponding metadata variable. Note that because the enrollment was ongoing throughout  the study period, and the treatment was taking more than a year, there are progressively fewer samples at the higher visit numbers. There were no study drop-outs.

The input driver script [examples/example_project_01/example_project.r](mgsat/src/master/examples/example_project_01/example_project.r) that the hypothetical user has created to run the analysis is included into the code repository along with the required Mothur count files and clinical metadata table. Comments within the script also describe how to install required R dependencies. You should also edit the location of the directory into which you checked out the MGSAT code. The analysis should run both on Linux and Windows. SNOW package is used for parallelization, with a default number of processes set to four.
The driver script expects to find its input files in the current working directory, and will generate the output report files in the same directory. The top-level HTML report file is called 0-report.html. This file can be opened directly from disk (use Firefox or Chrome; Internet Explorer has compatibility issues). The same directory also contains Markdown files that served as source of the HTML files. They can be converted to other formats such as Microsoft Word by using [Pandoc](http://johnmacfarlane.net/pandoc/).

You can open the [pre-computed HTML report](http://andreyto.bitbucket.org/mgsat/example_project_01/0-report.html) generated from running the example driver script (use Firefox or Chrome; Internet Explorer has compatibility issues).
The header section briefly describes how to navigate the report and subreports linked from its main page. Note that the report includes links to input files, as well as to various intermediate datasets in case user will want to re-analyze data in some custom ways. Each sections also describes the parameters used to run a specific method as well as references to corresponding R packages.

The user has analyzed data in subsets, such as controls vs patients before diet or patients before and after the start of the diet regimen. MGSAT placed the results for each subset into a separate subreport, linked from the top level page.

The example also shows how custom user analysis code `extra.method.task` can be injected into the pipeline and executed along with the already available methods. 

##Author##
Andrey Tovchigrechko `<andreyto AT gmail.com>`

##License##
GPLv3. See also COPYING file that accompanies the source code.

![Abundance profile patients across visits with and without drug treatment before dieting](http://andreyto.bitbucket.org/mgsat/example_project_01/plots/105879431cf3.png "Abundance profile patients across visits with and without drug treatment before dieting")