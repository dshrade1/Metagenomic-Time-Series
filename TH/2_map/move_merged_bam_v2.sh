cd v2_bam/
DATA_DIRS=*
for i in $DATA_DIRS; do
echo "$i"
mv $i/all_v2_$i.bam ../merged_bam_v2/all_v2_$i.bam
done

