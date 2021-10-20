#!/bin/sh

expList=("CanaryData/05232017_WT129_-6-7_NoQS")


for expPath in ${expList[@]};

do
echo $expPath

sh ./CanaryFileSort.sh "$expPath"

sh ./CanaryBkgndCorr.sh "$expPath"

sh ./MakeCanaryCorrectedMovie_remove_t1.sh "$expPath"

sh ./GetColonyPeaks.sh "$expPath"

sh ./CanaryColonyProcessing.sh "$expPath"

done