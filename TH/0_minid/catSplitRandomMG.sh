# catSplitRandomMG.sh TroutBogReads

cd temp

# concatentate random subsampled files
cat * > random_MG.fasta

# split into <20 MB files
mkdir random
split -l 120000 random_MG.fasta random/random

# list split read files (for mapping later)

BEGIN=0.7
INTERVAL=0.01
END=1.0

for a in `seq $BEGIN $INTERVAL $END`; do
DATA_DIRS=random/
for i in $DATA_DIRS; do
DATA_FILES="$i"*
for k in $DATA_FILES; do
echo "${k:0:6} ${k:7:12}" "$a" >> randomReadFileListNumbers.txt
done
done
done

# transfer metagenome to CHTC
sshpass -p "Geraldine2012" scp -r random dgshrader@submit-5.chtc.wisc.edu:~/"$1"



# transfer file list to CHTC
sshpass -p "Geraldine2012" scp -r randomReadFileListNumbers.txt dgshrader@submit-5.chtc.wisc.edu:~/"$1"
