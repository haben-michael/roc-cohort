#!/bin/bash

# use the current directory
#$ -cwd
#$ -S /bin/bash

# mail this address
# -M haben.michael@stanford.edu
# send mail on begin, end, suspend
# -m bes

/afs/ir/users/h/a/habnice/Desktop/R/R-3.2.1/bin/Rscript $*


## submit command
# for b in {1..500}
# do
#	qsub submit.sh $*
# done
