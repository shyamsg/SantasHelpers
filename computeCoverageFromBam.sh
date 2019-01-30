#! /bin/bash


if [ "x$1" == "x" ]; then 
    echo "Usage: computeCoverageFromBam.sh bamfilename"
    exit 1
fi
if [ ! -e $1 ]; then 
    echo "File $1 does not exist."
    exit 2
fi

inbam=$1
base=$(basename $1 .L.Dalen_14_wolf.scf.realigned.bam)

samtools stats -c 1,2000,1 $inbam > $base.stats
