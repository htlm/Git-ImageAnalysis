#!/bin/bash

expPath=$1
cd $expPath
echo "Making corrected Tiff files..."

mkdir -v rgb-corrected-noT0-tiff red-corrected-noT0-tiff green-corrected-noT0-tiff rgb-corrected-noT0-movie
cd ../..

/Applications/MATLAB_R2017a.app/bin/matlab -nodisplay -nodesktop -nosplash -r \
	"fnCanaryRGB_remove_t0001('$expPath');, exit"