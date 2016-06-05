#/usr/bin/time --verbose sh -c 'cat TH_fastq_ID_test2.txt | parallel -j 12 --workdir $PWD ./formatSubsetTransfer_fastq.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/ {}'


i="$2"

echo "$i"

# move to temp folder
cd temp

# copy read file to temp folder
cp "$1"/$i $i

# convert fastq to fasta in temp folder
python ../fastq2fasta.py $i
rm $i

# unwrap fasta
python ../FastaMLtoSL.py $i.fasta
mv $i.fasta.out $i.fasta

# length filter
perl ../removesmalls.pl 150 $i.fasta > $i.fasta.out

# subsample
python ../subsampler.py $i.fasta.out 0.02 > sub_$i.fasta

# remove the original read files and intermediate files from temp
rm $i.fasta.out
rm $i.fasta

# transfer to CHTC
# sshpass -p "Geraldine2012" scp -r */ dgshrader@submit-5.chtc.wisc.edu:~/"$2"
