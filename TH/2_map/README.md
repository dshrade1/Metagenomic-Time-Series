In this step, you map reads from individual metagenomes back to a combined assembly.

Input:  
metagenome read files  
coassembly  




After mapping:  
  
Move all the output bam files from their original folders to v2_bam, a folder of only bam files (1 folder of only all bam files per metagenome):  
Run the following from TroutBogReads in CHTC:  
./move_bam_v2.sh  
  
Merge the bam files into 1 bam file per metagenome.  
Run the following from TroutBogReads in CHTC:  
./merged_bam_v2.sh  
  
Move the merged bam files to a single folder:  
Run from TroutBogReads in CHTC:  
./move_merged_bam_v2.sh
