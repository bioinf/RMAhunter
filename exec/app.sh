#!/bin/bash
key=$1
coding=$2
maxafs=$3

./exec/hunter.py \
  -f /tmp/"$key".xvcf \
  -b /tmp/"$key".xbed \
  -c "$coding" \
  -m "$maxafs" \
  -o web/results/"$key" \
  -v False

cd web/results/$key/data
echo $(
  for i in 1 2 3 4; do
    split -l 30 -a4 -d tbl.0.t$i tbl.0.t$i.p
    [[ -f tbl.0.t$i.p0000 ]] || touch tbl.0.t$i.p0000
    cat tbl.0.t$i | wc -l
  done
)

# find ./E* -mtime +1 -delete
# find ./public/res/E* -mtime -1 -delete
