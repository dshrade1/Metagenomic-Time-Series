cd temp
while read i; do
echo "Receiving bam file" $i "..."
sshpass -p "Geraldine2012" scp dgshrader@submit-5.chtc.wisc.edu:~/TroutBogReads/merged_bam_gz/all_$i.bam.gz /home/dgshrader/Metagenomic-Time-Series/TH/3_count/temp
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
/home/dgshrader/Metagenomic-Time-Series/data/coassembly
# htseq-count -f sam -r name -s no -a 0 -t CDS -i locus_tag -m intersection-strict sorted_$i.sam ../../../data/coassembly/3300000553.gff > $i.count
htseq-count -f sam -r name -s no -a 0 -t repeat_region -i locus_tag -m intersection-strict sorted_$i.sam ../../../data/coassembly/3300000553.gff > $i_rr.count
# htseq-count -f sam -r name -s no -a 0 -t tRNA -i locus_tag -m intersection-strict sorted_$i.sam ../../../data/coassembly/3300000553.gff > $i_tRNA.count
# htseq-count -f sam -r name -s no -a 0 -t rRNA -i locus_tag -m intersection-strict sorted_$i.sam ../../../data/coassembly/3300000553.gff > $i_rRNA.count
rm sorted_$i.sam
rm sorted_$i.bam
rm reformatted_$i.bam
rm all_$i.bam
done < ../"$1"

#nohup ./receive_gz.sh readFileList_aa &
#nohup ./receive_gz.sh readFileList_ab &
#nohup ./receive_gz.sh readFileList_ac &
#nohup ./receive_gz.sh readFileList_ad &
#nohup ./receive_gz.sh readFileList_ae &
#nohup ./receive_gz.sh readFileList_af &
#nohup ./receive_gz.sh readFileList_ag &
#nohup ./receive_gz.sh readFileList_ah &
