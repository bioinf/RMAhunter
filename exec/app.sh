#!/bin/bash
key=$1
coding=$2
maxafs=$3
./build/exec/hunter /tmp/$key.xvcf /tmp/$key.xbed \
 "$coding" "$maxafs" build/web/results/$key \
 build/data/sdf.csv build/data/sdf_plus.csv

cd build/web/results/
echo $(
  for i in 1 2 3 4; do
    cat $key.t$i | sort | awk {'print $2'} > $key.ts$i
    split -l 30 -a3 -d $key.ts$i $key.ts$i.p
    [[ -f $key.ts$i.p000 ]] || touch $key.ts$i.p000
    cat $key.ts$i | wc -l
  done
)

# find ./E* -mtime +1 -delete
# find ./public/res/E* -mtime -1 -delete
