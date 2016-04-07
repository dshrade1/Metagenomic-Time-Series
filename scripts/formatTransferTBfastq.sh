cd "$1"
DATA_PATH_FILES=*.fastq
for i in $DATA_PATH_FILES; do
echo "$i"

# copy read file to temp folder
cp $i /home/dgshrader/Metagenomic-Time-Series/data/temp/$i


# convert fastq to fasta in temp folder
cd /home/dgshrader/Metagenomic-Time-Series/data/temp/
# awk 'BEGIN{P=1}{if(P==1||P==2){gsub(/^[@]/,">");print}; if(P==4)P=0; P++}' $i > $i.fasta
python /home/dgshrader/Metagenomic-Time-Series/scripts/fastq2fasta.py $i

# unwrap fasta
python ~/Metagenomic-Time-Series/scripts/FastaMLtoSL.py $i.fasta
mv $i.fasta.out $i.fasta

# length filter
perl ~/Metagenomic-Time-Series/scripts/removesmalls.pl 150 $i.fasta > $i.fasta.out

# split
elevenletter=${i:0:11}
mkdir /home/dgshrader/Metagenomic-Time-Series/data/temp/${elevenletter}
split -l 120000 /home/dgshrader/Metagenomic-Time-Series/data/temp/$i.fasta.out /home/dgshrader/Metagenomic-Time-Series/data/temp/${elevenletter}/${elevenletter}

# transfer to CHTC
sshpass -p "Geraldine2012" scp -r /home/dgshrader/Metagenomic-Time-Series/data/temp/*/ dgshrader@submit-5.chtc.wisc.edu:~/TroutBogReads
rm -r /home/dgshrader/Metagenomic-Time-Series/data/temp/*/
rm *.fasta.out
rm *.fastq
rm *.fastq.fasta
cd "$1"
done

