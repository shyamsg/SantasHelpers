#! /bin/bash

isComputerome=$(module avail angsd 2>&1 | grep -c angsd)

if [ "x$isComputerome" == "x0" ]; then 
    ANGSD="$HOME/.local/bin/angsd"
    REALSFS="$HOME/.local/bin/realSFS"
else
    module load angsd
    ANGSD=$(which angsd)
    REALSFS=$(which realSFS)
fi

inbam=$1
reffasta=$2
bootstrap=$3
threads=$4
extn=$(basename $inbam | cut -f2- -d.)
out=$(basename $inbam .$extn)

if [ ! -s $out.saf.gz ]; then 
    $ANGSD -i $inbam -anc $reffasta -dosaf 1 -fold 1 -gl 1 -out $out -nThreads $threads
fi
#followed by the actual estimation
if [ ! -s $out.sfs ]; then
    $REALSFS $out.saf.idx -bootstrap $bootstrap -P $threads > $out.sfs
fi
