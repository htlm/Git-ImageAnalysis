#!/bin/bash

echo $1
expPath=$1


/Applications/MATLAB_R2018a.app/bin/matlab -nodisplay -nodesktop -nosplash -r \
	"tic, fnCanaryColonyProcessingNoT0('$expPath', [1 288], 1); toc , exit"