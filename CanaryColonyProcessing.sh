#!/bin/bash

echo $1
expPath=$1

/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -nodesktop -nosplash -r \
	"fnCanaryColonyProcessingNoT0('$expPath', [1, 288] , 1);, exit"