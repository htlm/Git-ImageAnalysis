# Git-ImageAnalysis

This code requires MATLAB2018a
The default path to the MATLAB installation is '/Applications/MATLAB_R2018a.app/bin/matlab'. 
If this is not the correct path, the path will need to be adjusted in the following files to use the pipeline:
- CanaryFileSort.sh 
- CanaryBkgndCorr.sh 
- MakeCanaryCorrectedMovie_remove_t1.sh 
- GetColonyPeaks.sh 
- CanaryColonyProcessing.sh

To run this code on the provided dataset, run the file 'CanaryProcessingFullPipeline.sh' from a bash shell. 

After a successful run of the pipeline, the processed data requires an additional 750MB of memory. 
