#! /bin/bash


##############################################################
# PBS system usage of admixture to run replicates as a single#
# array job. Note that the usage is (all 3 args mandatory)   #
# runAdmixtureArray.sh inbed K numThreads                    #
##############################################################

module load admixture/1.23

rep=$PBS_ARRAYID
mkdir -p rep$rep

K=$2
numThreads=$3
inbed=$1
inbim=${inbed/.bed/.bim}
infam=${inbed/.bed/.fam}

ln -s $inbed .
ln -s $inbim .
ln -s $infam .

infile=`basename $inbed`

admixture -j$numThreads --seed=$RANDOM --cv $infile $K 
