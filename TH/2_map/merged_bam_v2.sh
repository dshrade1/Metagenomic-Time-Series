cd v2_bam
DATA_DIRS=*
for i in $DATA_DIRS; do
echo "$i"
samtools merge $i/all_v2_$i.bam $i/*.bam
done

