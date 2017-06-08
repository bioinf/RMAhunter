#!/bin/bash
key=$1
gsx=$2
newkey=$3

cd web/results/
mkdir $newkey
mkdir $newkey/data

echo $(
  for i in 1 2 3 4; do
    cat $key/data/tbl.0.t$i | grep -P $gsx > $newkey/data/tbl.0.t$i
    split -l 30 -a4 -d $newkey/data/tbl.0.t$i $newkey/data/tbl.0.t$i.p
    [[ -f $newkey/data/tbl.0.t$i.p0000 ]] || touch $newkey/data/tbl.0.t$i.p0000
    cat $newkey/data/tbl.0.t$i | wc -l
  done
)
