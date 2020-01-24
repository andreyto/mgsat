#/bin/sh
set -ex
sample_id_rx="$1"
[ -n "$sample_id_rx" ] || exit 1
shift
sample_dir="$1"
[ -n "$sample_dir" ] || exit 1


this_dir=$(dirname $0)
dada_dir=$(cd "$this_dir"/.. && pwd)
wdl_dir=$(cd "$dada_dir"/wdl && pwd)

sample_sheet_file=sample_sheet.csv

threads=32

## if given a directory name, try to turn it into abs path,
## but it can also be a wildcard
if [ -d "$sample_dir" ]; then
  sample_dir=$(cd "$sample_dir" && pwd)
fi

$dada_dir/make_sample_sheet.sh "$sample_id_rx" "$sample_sheet_file" "$sample_dir/*.gz"
$wdl_dir/run_dada_wdl.sh

