#!/bin/bash
key=$1
gsx=$2
newkey=$3

cd ./public/res/

echo $(
  for i in 1 2 3 4; do
    cat $key.t$i | grep -P $gsx | sort | awk {'print $2'} > $newkey.ts$i
    split -l 30 -a3 -d $newkey.ts$i $newkey.ts$i.p
    [[ -f $newkey.ts$i.p000 ]] || touch $newkey.ts$i.p000
    cat $newkey.ts$i | wc -l
  done
)
