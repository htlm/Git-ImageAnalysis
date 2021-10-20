expPath=$1
cd $expPath
echo "File Sorting..."

mkdir -v DsRed-DsRed GFP-GFP DsRed-none GFP-none red-corrected green-corrected
mv *em_DsRed* DsRed-DsRed
mv *em_GFP* GFP-GFP
mv *DsRed_em_None* DsRed-None
mv *GFP_em_None* GFP-None

cd ../..
echo "Done"