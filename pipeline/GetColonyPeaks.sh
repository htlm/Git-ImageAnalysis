#!/bin/bash

echo $1
expPath=$1

/Applications/MATLAB_R2018a.app/bin/matlab -nodisplay -nodesktop -nosplash -r \
	"fnPkFndSurvey('$expPath'), exit"