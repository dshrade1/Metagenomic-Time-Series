Run all the following from ~/Metagenomic-Time-Series/TH/0_minid: 


1. To import, format, and subsample the FASTA files:
/usr/bin/time --verbose sh -c 'cat TH_fasta_ID.txt | parallel -j 12 --workdir $PWD ./formatSubsetTransfer_fasta.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/ {}'

2. To import, format, and subsample the FASTQ files:
/usr/bin/time --verbose sh -c 'cat TH_fastq_ID.txt | parallel -j 15 --workdir $PWD ./formatSubsetTransfer_fastq.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/fastq {}'

3. to concatenate the subsampled files into a single file and send to CHTC (with a randomReadFileList.txt):
./catSplitRandomMG.sh TroutBogReads

