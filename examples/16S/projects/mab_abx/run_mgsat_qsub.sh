#!/bin/sh
echo "Rscript ~/work/mgsat/examples/16S/projects/mab_abx/mgsat_project.R" | qsub -V -d `pwd` -S /bin/sh -l nodes=1:ppn=16

