# time ./formatTransferTBfasta_test.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/ TroutBogReads TH_fasta_ID_test.txt

while read i; do
echo "$i"

# move to folder containing the data
#cd "$1"

# move to temp folder
cd temp

# copy read file to temp folder
cp "$1"/$i $i

# unwrap fasta
python ../FastaMLtoSL.py $i
# it puts the i.out file in the temp folder, where it got its i file.
rm $i
mv $i.out $i

# split read file
fourletter=${i:0:4}
mkdir ${fourletter}
split -l 140000 $i ${fourletter}/${fourletter}

# list the new files that were created (for mapping later)
DATA_DIRS=*/
for i in $DATA_DIRS; do
DATA_FILES="$i"*
for k in $DATA_FILES; do
echo "${k:0:4} ${k:5:11}" >> readFileList.txt
done
done

# transfer read folder to CHTC
sshpass -p "Geraldine2012" scp -r /home/dgshrader/Metagenomic-Time-Series/TH/1_preprocess/temp/*/ dgshrader@submit-5.chtc.wisc.edu:~/"$2"

# remove read files from temp
rm -r /home/dgshrader/Metagenomic-Time-Series/TH/1_preprocess/temp/*/
rm /home/dgshrader/Metagenomic-Time-Series/TH/1_preprocess/temp/*.fasta
done < "$3"

# consider: see if "$1"/$i works or "$1"/"$i"
