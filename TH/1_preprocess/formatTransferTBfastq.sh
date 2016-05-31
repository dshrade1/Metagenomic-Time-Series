# time ./formatTransferTBfastq_test.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/fastq/ TroutBogReads TH_fastq_ID_test.txt

while read i; do
echo "$i"

# move to folder containing the data
#cd "$1"

# move to temp folder
cd temp

# copy read file to temp folder
cp "$1"/$i $i

# convert fastq to fasta in temp folder
# awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $i > $i.fasta
python ../fastq2fasta.py $i

# unwrap fasta
python ../FastaMLtoSL.py $i.fasta
mv $i.fasta.out $i.fasta

# length filter
perl ../removesmalls.pl 150 $i.fasta > $i.fasta.out

# split
elevenletter=${i:0:11}
mkdir ${elevenletter}
split -l 120000 $i.fasta.out ${elevenletter}/${elevenletter}

# list split read files (for mapping later)
DATA_DIRS=*/
for i in $DATA_DIRS; do
DATA_FILES="$i"*
for k in $DATA_FILES; do
echo "${k:0:11} ${k:12:23}" >> readFileList.txt
done
done

# transfer to CHTC
sshpass -p "Geraldine2012" scp -r */ dgshrader@submit-5.chtc.wisc.edu:~/"$2"
rm -r */
rm *.fasta.out
rm *.fastq
rm *.fastq.fasta
done < "$3"

sshpass -p "Geraldine2012" scp -r readFileList.txt dgshrader@submit-5.chtc.wisc.edu:~/"$2"
#rm readFileList.txt 
