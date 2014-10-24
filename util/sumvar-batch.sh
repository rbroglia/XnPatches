#!/bin/bash -
set -o nounset                              # Treat unset variables as an error

ZNUMS=$(head -n 1 sumvar-batch.par | awk '{print $1}')
VNUM=$(tail -n 1 sumvar-batch.par | awk '{print $1}')
FOUTNAME="sum-var-$VNUM.dat"

sed s/SEDZNUMS/$ZNUMS/ sumvar-batch.mcr |  sed s/SEDVNUM/$VNUM/ |  sed s/SEDFOUTNAME/$FOUTNAME/ > tmp.mcr
