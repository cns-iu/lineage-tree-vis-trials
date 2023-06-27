#!/bin/bash
source constants.sh
set -ev

for measure in betweenness pagerank; do
  for dsName in $(ls ./input-data)
  do
    network=./input-data/$dsName
    dataset=${dsName}_$measure
    echo $network $dataset
    echo
    OUT=./datasets/${dataset}_BatchTree
    OUT2=./datasets/${dataset}_DELG
    OUT3=./datasets/${dataset}_CG
    mkdir -p $OUT $OUT2 $OUT3

    if [ ! -e $OUT/network.dot ]
    then
      MAPPING=$network/${dsName}_labels.csv
      NW_FILE=$network/${dsName}_tree.nw
      python3 src/newick2dot.py --mapping $MAPPING --node-weights $measure \
        $NW_FILE $OUT/network.dot &> $OUT/convert.log.txt
      cp $OUT/network.dot $OUT2/network.dot
      cp $OUT/network.dot $OUT3/network.dot
    fi
    if [ ! -e $OUT/config.sh ]
    then
      cat src/config-template.sh | perl -pe "s/--DATASET--/${dataset}_BatchTree/g" > $OUT/config.sh
      cat src/config-template.sh | perl -pe "s/--DATASET--/${dataset}_DELG/g;s/BatchTree/DELG/g" > $OUT2/config.sh
      cat src/config-template.sh | perl -pe "s/--DATASET--/${dataset}_CG/g;s/BatchTree/CG/g" > $OUT3/config.sh
    fi
  done
done
