#!/bin/bash
source constants.sh
set -ev

for network in $(ls ./input-data/*.txt)
do
  dataset=$(basename "${network%.txt}")
  echo $network $dataset
  echo
  OUT=./datasets/${dataset}_BatchTree
  OUT2=./datasets/${dataset}_DELG
  OUT3=./datasets/${dataset}_CG
  mkdir -p $OUT $OUT2 $OUT3

  if [ ! -e $OUT/network.dot ]
  then
    python3 src/newick2dot.py ${network} $OUT/network.dot &> $OUT/convert.log.txt
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
