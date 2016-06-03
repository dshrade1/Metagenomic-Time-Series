DATA_DIRS=*/
for i in $DATA_DIRS; do
echo "$i"
mkdir v2_bam/$i
mv $i/v2_* v2_bam/$i/
done
