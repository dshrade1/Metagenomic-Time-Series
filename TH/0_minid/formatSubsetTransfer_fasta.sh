# /usr/bin/time --verbose sh -c 'cat TH_fasta_ID.txt | parallel -j 12 --workdir $PWD ./formatSubsetTransfer_fasta.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/ {}'

i=$2

echo "$i"

# move to temp folder
cd temp

# copy read file to temp folder
cp "$1"/$i $i

# unwrap fasta
python ../FastaMLtoSL.py $i
rm $i
mv $i.out $i

# subsample
python ../subsampler.py $i 0.02 > sub_$i

# remove the original read files from temp
rm $i




# transfer read folder to CHTC
# sshpass -p "Geraldine2012" scp -r /home/dgshrader/Metagenomic-Time-Series/TH/1_preprocess/temp/*/ dgshrader@submit-5.chtc.wisc.edu:~/"$2"
