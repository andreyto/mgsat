#/bin/sh
[ -n "$1" ] || exit 1
[ -n "$2" ] || exit 1
[ -n "$3" ] || exit 1
python -m MICGENT.workflow_util iterate-sequence-run \
    --samp-id-extractor "$1" \
    --out-file "$2" \
    --out-file-format csv "$3" \
    --samp-id-nomatch warn

