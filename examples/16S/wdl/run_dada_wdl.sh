#!/bin/sh
set -ex
this_dir=$(cd $(dirname $0) && pwd)
export PATH=$this_dir/..:$PATH
rm -f wdls.zip; zip -j wdls.zip $this_dir/*.wdl
cp $this_dir/../bbduk_adapters.fa ./
cromwell \
    -Xms64g -Xmx64g \
    -Dconfig.file=$this_dir/wdl.conf \
    run \
    --inputs $this_dir/dada_inputs.json \
    --options $this_dir/mtg16s_wf_opt.json \
    --imports wdls.zip \
    $this_dir/dada.wdl  \
    2>cromwell.err 1>cromwell.out
#find wdl_final/ -name '*.fasta.gz' | grep 'collect_outputs/execution/glob'

