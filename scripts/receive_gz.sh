cd /home/dgshrader/metagenomic_ts/data/temp
while read i; do
echo "Receiving bam file" $i "..."
sshpass -p "Geraldine2012" scp dgshrader@submit-5.chtc.wisc.edu:~/TroutBogReads/merged_bam_gz/all_$i.bam.gz /home/dgshrader/metagenomic_ts/data/temp
echo "Opening gz file..."
gunzip all_$i.bam.gz
# only run the following two lines if the cigar strings are not properly formatted in the bam file (i.e., if you did not include "sam=1.3" in run_bbmap.sh
# if the cigar strings are improperly formatted, also use the hashed samtools sort line and comment out the currently active samtools sort line.
#echo "Reformatting bam file" $i "..."
#/home/dgshrader/metagenomic_ts/scripts/bbmap/reformat.sh in=all_$i.bam out=reformatted_$i.bam sam=1.3 -Xmx6g t=8
echo "Sorting bam file" $i "..."
#samtools sort -l 0 -n -m 6G reformatted_$i.bam sorted_$i
samtools sort -l 0 -n -m 6G all_$i.bam sorted_$i
echo "Converting to sam..."
samtools view -h sorted_$i.bam > sorted_$i.sam
echo "Counting reads in" $i "..."
htseq-count -f sam -r name -s no -a 0 -t CDS -i locus_tag -m intersection-strict sorted_$i.sam /home/dgshrader/metagenomic_ts/data/coassembly/3300000553.gff > $i.count
rm sorted_$i.sam
rm sorted_$i.bam
rm reformatted_$i.bam
rm all_$i.bam
done < /home/dgshrader/metagenomic_ts/data/metagenomic_reads/readFileList_test.txt
#nohup ./receive_gz.sh &
