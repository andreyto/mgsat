# 16S analysis files for the submitted manuscript

**Jones et al, Antibacterial Monoclonal Antibodies Do Not Disrupt the Intestinal Microbiome or its Function**

- `run_data.sh` - executes primary analysis making taxonomic abundance 
  matrix from FASTQ files
- `run_mgsat_qsub.sh` - executes secondary analysis starting from the
  taxonomic abundance matrix and study design table

Run these in separate directories for each of the three experiments.
For the primary analysis, FASTQ files can be downloaded from the ENA: 
https://www.ebi.ac.uk/ena/data/view/PRJEB34462
You need the Submitted Files flavor that retains the original file names
because the sample IDs will be parsed out of the file names.

The pre-computed outputs of the primary analysis are saved in the subdiretories
here, so that you could run the secondary analysis without redoing the primary
one.

Note that the manuscript relabelled some treatment groups for clarity and 
consistency as compared to the labels used in the StudyDesign.txt files here
and in the ENA file names:

| From  | To     |
| ----- | ------ |
| R347  | c-IgG  |
| Naive | Saline |
| LC10  | 4893   |

The secondary analysis will create the graphical reports as described
in the README at the root of this repository.

