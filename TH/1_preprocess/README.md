Purpose:  
Format metagenomic reads for mapping, send over to CHTC.  
  
Input:  
TH_fasta_ID.txt # list of metagenome read fasta files to map in CHTC  
TH_fastq_ID.txt # list of metagenome read fastq files to map in CHTC  
/data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog # directory containing TB metagenome fasta files  
/data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/fastq # directory containing TB metagenome fastq files  
  
Scripts:  
formatTransferTBfasta.sh # formats fastq metagenome read files and sends folder of processed files to CHTC  
Note: this script calls FastaMLtoSL.py  
formatTransferTBfastq.sh # formats fasta metagenome read files and sends folder of processed files to CHTC  
# Note: this script calls fastq2fasta.py, FastaMLtoSL.py, and removesmalls.pl  
  
Software:  
Python 2.7  
Perl v5.18.2  
  
Result:  
One folder per metagenome is created and sent over to CHTC.  
  
Usage:  
  
./formatTransferTBfasta.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog TH_fasta_ID.txt  
./formatTransferTBfastq.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/fastq TH_fastq_ID.txt  
  
Memory: TBD  
  
Disk: TBD  
  
Time: TBD  
