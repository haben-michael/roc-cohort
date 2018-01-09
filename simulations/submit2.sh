#!/bin/bash
#SBATCH --time=40
#SBATCH --nodes=1
#SBATCH --mem=1GB
#SBATCH -c 1
# use the current directory
#$ -cwd
#$ -S /bin/bash

# mail this address
# -M haben.michael@stanford.edu
# send mail on begin, end, suspend
# -m bes

/home/habnice/R/bin/Rscript $*


## submit command
# for b in {1..500}
# do
#	qsub submit.sh $*
# done
