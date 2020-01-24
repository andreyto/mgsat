#!/bin/sh
set -ex
cromwell -Dconfig.file=wdl.conf run ~/work/mgsat/examples/16S/wdl/bbduk.wdl  bbduk_inputs.json bbduk_wf_opt.json
find wdl_final/ -name '*.fasta.gz' | grep 'collect_outputs/execution/glob'

