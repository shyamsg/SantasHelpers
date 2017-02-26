#! /bin/bash 


numRegions=$1
regionSize=$2
outName=$3
inbamList=$4
infai=$5


### Make the regions file from bedtools or shyam script
bedtools random -n $numRegions -l $regionSize -g $infai | awk '{ print $1:($2-1)"-"($3-1);}'> $outName.chosen.bed

### Make the consensus sequences, run trimal on it and them run fasttree
echo "#! /bin/bash

module unload samtools
module load samtools/1.2.1
module load trimal/1.4.1




"


