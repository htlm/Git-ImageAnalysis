#!/bin/sh

echo $1
expPath=$1

/Applications/MATLAB_R2018a.app/bin/matlab -nodisplay -nodesktop -nosplash -r \
	"fnDatasetBkgndCorrectionUpdated('$expPath' , 1:288 );, exit"