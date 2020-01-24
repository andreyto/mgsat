#/bin/bash
set -ex
mock="$1"

this_dir=$(dirname $0)
dada_dir=$(cd "$this_dir"/../.. && pwd)
wdl_dir=$(cd "$dada_dir"/wdl && pwd)

sample_dir=/path/to/FASTQ

if [ -n "$mock" ]; then
    sample_id_rx='^(M[^_]+_[^_]+)_'
else
    sample_id_rx='^(D[^_]+)_'
fi
$wdl_dir/run_dada_project.sh "$sample_id_rx" "$sample_dir"

