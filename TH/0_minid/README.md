Run all the following from ~/Metagenomic-Time-Series/TH/0_minid: 


1. To import, format, and subsample the FASTA files:
/usr/bin/time --verbose sh -c 'cat TH_fasta_ID.txt | parallel -j 12 --workdir $PWD ./formatSubsetTransfer_fasta.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/ {}'

2. To import, format, and subsample the FASTQ files:
/usr/bin/time --verbose sh -c 'cat TH_fastq_ID.txt | parallel -j 15 --workdir $PWD ./formatSubsetTransfer_fastq.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/fastq {}'

3. to concatenate the subsampled files into a single file and send to CHTC (with a randomReadFileListNumbers.txt):
./catSplitRandomMG.sh TroutBogReads

4. in CHTC, run the mapping.
condor_submit run_bbmap_random.sub
Note: this calls run_bbmap_random.sh

5. After mapping, in CHTC run listReadFiles_random.sh to generate a list readFileList_random.txt, which is an input into percent_mapped.sh.

6. Still in CHTC, run percent_mapped.sh to generate a file of percentages mapped per minid value tested. Note: I could also investigate

7. Transfer files from CHTC to local machine for analysis & plotting in R.

8. Run select_minid.R to produce plots of the % reads mapped and % bases mapped across a range of different values of minid
