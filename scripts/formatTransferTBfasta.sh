# time ./formatSplitTransfer.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/ 
cd "$1"
DATA_PATH_FILES=*.fasta
for i in $DATA_PATH_FILES; do
echo "$i"

# copy read file to temp folder
cp $i /home/dgshrader/Metagenomic-Time-Series/data/temp/$i

# unwrap fasta
python /home/dgshrader/Metagenomic-Time-Series/scripts/FastaMLtoSL.py /home/dgshrader/Metagenomic-Time-Series/data/temp/$i
# it puts the i.out file in the temp folder, where it got its i file.
rm /home/dgshrader/Metagenomic-Time-Series/data/temp/$i
mv /home/dgshrader/Metagenomic-Time-Series/data/temp/$i.out /home/dgshrader/Metagenomic-Time-Series/data/temp/$i

# split read file
fourletter=${i:0:4}
mkdir /home/dgshrader/Metagenomic-Time-Series/data/temp/${fourletter}
split -l 140000 /home/dgshrader/Metagenomic-Time-Series/data/temp/$i /home/dgshrader/Metagenomic-Time-Series/data/temp/${fourletter}/${fourletter}

# transfer to CHTC
sshpass -p "Geraldine2012" scp -r /home/dgshrader/Metagenomic-Time-Series/data/temp/*/ dgshrader@submit-5.chtc.wisc.edu:~/TroutBogReads
rm -r /home/dgshrader/Metagenomic-Time-Series/data/temp/*/
rm /home/dgshrader/Metagenomic-Time-Series/data/temp/*.fasta
done
