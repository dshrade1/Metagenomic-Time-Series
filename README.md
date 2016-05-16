## Introduction

This is a protocol to map metagenomic reads to a metagenomic coassembly. Steps covered are as follows:  

1. Set up directory for pre-mapping formatting in lab server (Zissou)  
2. Set up directory for mapping in compute cluster (CHTC)  
3. Transfer coassembly file to CHTC  
4. Format and transfer metagenomic read files to CHTC  
5. Prep CHTC submit node for mapping  
6. Test mapping in CHTC  
7. Run mapping in CHTC  
8. Verify mapping.  
9. Merge mapping output by metagenome.  
10. Compress mapping output for transfer to Zissou.  
11. Import mapping output to Zissou; count reads per gene.  
12. Quantify reads mapped per COG annotation.  
  
##**Step 1: Set up directory in Zissou for pre-mapping formatting**##

**1a. Pull this repo into your Zissou folder or otherwise set up your Zissou folder the way this repo is structured:**  
  
Metagenomic-Time-Series  
|-- scripts  
|--|-- FastaMLtoSL.py             # called by formatTransfer scripts  
|--|-- fastq2fasta.py             # called by formatTransferTBfastq.sh  
|--|-- formatTransferTBfasta.sh   # formats TB reads (fasta), sends to CHTC  
|--|-- formatTransferTBfastq.sh   # formats TB reads (fastq), sends to CHTC  
|--|-- removesmalls.pl            # called by formatTransferTBfastq.sh  
|-- scripts_CHTC  
|--|-- run_bbmap.sh               # executable for mapping in CHTC  
|--|-- run_bbmap.sub              # submit file for mapping in CHTC  
|--|-- listReadFiles1.sh          # creates list of files to map  
|--|-- listReadFiles2.sh          # creates list of files to map  
|-- data  
|--|-- coassembly  
|--|-- temp  
|-- README.md  

**1b. Make sure the following are installed in Zissou:**
- Python (I used 2.7)
- Perl (I used v5.18.2)
- sshpass (I used v1)

##**Step 2: Set up directory for mapping in CHTC**

**2a. Obtain CHTC account by filling out http://aci.wisc.edu/large-scale-request/.** You may meet with a CHTC rep (such as Christina Koch or Lauren Michael) who will provide you with an account. You'll need to get a SQUID account as well. **Note:** You may need to request an increase in your submit node folder's disk space quota to account for the total expected size of input and output files. I got 500 GB space total, which was just right for my purposes. Also, CHTC restricts transfer of files >20MB between submit nodes and compute nodes, so metagenomic read files will be split into <20MB.

**2b. Download CyberDuck or another SFTP tool to transfer files among your server (or local machine) and CHTC.** Set up the connection with your CHTC home folder selecting “SFT (SSH File Transfer Protocol)” rather than the default “FTP (File Transfer Protocol)”. Connection should look something like dgshrader@submit-5.chtc.wisc.edu. Set up your Zissou home directory in CyberDuck as well. 

**2c. In CHTC submit node, create a folder for the mapping project.** Here, let's call it **TroutBogReads**. Set it up as follows:  

TroutBogReads  
|-- bbmap_only   # this is a subset the bbmap software. This was from BBMap_35.82.tar.gz; updates have been made since.  
|--|-- bbmap.sh  
|--|-- build.xml  
|--|-- calcmem.sh  
|--|-- current    # this is the entire bbmap 'current' folder.  
|--|-- jni        # this is the entire bbmap 'jni' folder.  
|-- run_bbmap.sh  
|-- run_bbmap.sub  
|-- TB_Epi_ID.txt  
|-- listReadFiles1.sh  
|-- listReadFiles2.sh  

##**Step 3. Prep CHTC SQUID folder for mapping.** 
Add the .fna file for the combined assembly (hereafter, "coassembly") from Zissou to your squid folder on CHTC. 

Coassembly locations on Zissou are as follows:  
```
/data_lakes/Metagenomes/TroutBog/coassembly/3300000553/3300000553.a.fna # Trout Bog Hypolimnion, 537M  
/data_lakes/Metagenomes/TroutBog/coassembly/3300000439/3300000439.a.fna # Trout Bog Epilimnion, 255M  
/data_lakes/Metagenomes/Mendota/coassembly/3300002835/3300002835.a.fna # Mendota (Epilimnion?), 4.2G  
```
To get coassembly from Zissou to CHTC, first compress the coassembly .fna file and save this compressed file in your personal folder on Zissou:  


```
gzip -c /data_lakes/Metagenomes/TroutBog/coassembly/3300000553/3300000553.a.fna > /home/dgshrader/Metagenomic-Time-Series/data/coassembly/THcoassembly_100percent.fna.gz #this took ~3min for TB Hypolimnion
```

You can then transfer the gzipped coassembly by either:  

**(a) using Cyberduck:** Drag THcoassembly_100percent.fna.gz from /home/dgshrader/Metagenomic-Time-Series/data/coassembly/ to CHTC SQUID folder.  
OR  
**(b) scp:**  

```
scp /home/dgshrader/Metagenomic-Time-Series/data/coassembly/THcoassembly_100percent.fna.gz dgshrader@submit-5.chtc.wisc.edu:~/../../squid/dgshrader  
```

##**Step 4. Format read files in Zissou & transfer them to CHTC**

There are three types of metagenomic read files in Zissou, and the required formatting is different for each type.  

**(a) fasta files in /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/**  
- These are merged reads.  
- The text is wrapped and needs unwrapping.  
- These files **have** been length-filtered to contain only reads >150bp.  
- These files are identified by a 4-letter code before the read file extension.  
- To format and transfer these TB fasta files, from ~/Metagenomic-Time-Series/scripts type:  
```
time ./formatTransferTBfasta.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/
```

**(b) fastq files in /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/fastq.**  
- These are merged reads.  
- The text is wrapped and needs unwrapping.  
- These files have **not** been length-filtered to contain only reads >150bp, and need filtering.  
- These files need to be converted to fasta.  
- These files are identified by an 11-letter code before the read file extension.  
- To format and transfer these TB fastq files, from ~/Metagenomic-Time-Series/scripts type:  

```
time ./formatTransferTBfastq.sh /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/fastq/
```

**(c) fasta files in /data_lakes/Metagenomes/Mendota/MergedReads-Mendota/fasta/lenfiltered150**
- These are merged reads.  
- The text is **not** wrapped and does not need unwrapping.  
- These files **have** been length-filtered to contain only reads >150bp, and need no filtering.
- These files are identified by a 4-letter code before the read file extension.  
- The formatting script for these ME fastq files is still in development. See description for working with Mendota coassembly mapping at the end of this document.  

**Note:** maybe you want to use this protocol to map a set of reads that are none of the above 3. Ask:
Are they length-filtered? To check, type:  

```
awk '{ print length($0); }' /data_lakes/Metagenomes/TroutBog/MergedReads-TroutBog/IHPI.fasta > character_count.txt  
```
Are they wrapped or unwrapped text? To check, type:  
```
sed -n 1,2p IHXP.fasta > test.txt  
```
This shows you that each line is considered a new line, rather than 2 lines per read, as in the Mendota fasta files.  
You can mix and match the appropriate portions of the formatTransferTBfasta.sh and formatTransferTBfastq.sh based on the formatting needs for your files. In the end what you want is an unwrapped, length-filtered fasta file that has been split into <20MB pieces stored in a folder with the same code name as you find before the file extension of the read file. Scripts also use different approaches based on whether the code name before the file extension is 4 letters or 11 letters long, so account for this when modifying scripts.  

##**Step 5: Prep CHTC submit node for mapping**

The following steps should be performed after the metagenomic read files have been transferred to CHTC. By now, your split-up read files should be in a set of folders in the CHTC TroutBogReads folder, with one folder per metagenome. However, ALL Trout Bog metagenomes are here - not just epilimnion or just hypolimnion. So if you're mapping hypolimnion reads to the hypoliminion coassembly, you'll want to get rid of the epilimnion metagenomic read folders and vice versa.  

**5a. Delete metagenomic read files from non-target lake layer.**  
```
xargs rm -rf <TB_Epi_ID.txt
```

**5b. Create text file listing each metagenomic read folder name and read file name separated by a space.** Start by temporarily moving bbmap_only out of the TroutBogReads folder into your main CHTC submit node user folder (e.g., /home/dgshrader/bbmap_only). Run the following from TroutBogReads:  
```
mv bbmap_only ..
```

To create the read file lists, run the following lines from within the TroutBogReads folder.  

```
./listReadFiles1.sh
./listReadFiles2.sh
```
This will produce 2 text documents (readFileList1.txt and readfilelist2.txt) that each have some wrong characters in them. Download (Cyberduck) the resulting files readFileList1.txt and readfilelist2.txt, copy the appropriately-written lines of each one into readFileList.txt, and bring readFileList.txt into the CHTC TroutBogReads folder.  

The result should look like this:  
```
H0586_93745 H0586_93745aa  
H0586_93745 H0586_93745ab  
H0586_93745 H0586_93745ac  
.  
.  
.  
IHXY IHXYek  
IHXY IHXYel  
IHXY IHXYem  
```
Bring bbmap_only back into the TroutBogReads folder on CHTC:  
```
mv bbmap_only/ TroutBogReads/  
```

##**Step 6: Test mapping in CHTC**  

**6a. Find the largest metagenomic read files on CHTC.**  

In CHTC, run maxfilesizes.sh. Then download the output maxfilesizes.txt using Cyberduck. Open with Excel. Sort. Select the largest file names.  



**6b. Copy the first two lines of readFileList.txt into a document called readFileList_test.txt.** It should look like this:  

```
H0596_93745 H0596_93745au
H0596_93770 H0596_93770an
```

**6c. Edit the script "run_bbmap.sub"** to have the last line say  

```
queue dir,file from readFileList_test.txt
```

**6d. Run the test mapping: **  
```
condor_submit run_bbmap.sh  
```


**6d. Check memory and disk usage.**  
Within each metagenome's respective folder, there are log and error files for each <20MB piece. Open the log file and scroll to the bottom. This displays your memory and disk usage. Make sure your usage is under the amount you requested! If not, increase the disk and memory limits by editing the request_memory and request_disk lines in run_bbmap.sub. Open the .err files. Make sure othing strange happened (What could happen? See here.)  If you're good, proceed. Note: requesting more memory than needed could make you have to wait longer in the Condor queue, so you want to request enough that your jobs don't fail, but not so much memory that it delays your jobs a bunch in the queue.

##**Step 7:  Run the whole mapping**

**7a. Edit the script "run_bbmap.sub"** to have the last line say  

```
queue dir,file from readFileList.txt
```

**7b. Run the whole mapping**  
```
condor_submit run_bbmap.sh  
```

So you'll be running several mappings in parallel. To test the progress of the mapping, you can type:
```
condor_q $USER
```
This will list the remaining mappings and their progress.  
  
Congratulations! You've mapped the metagenomic reads to the coassembly.

##**Step 8: Check that the mappings worked.**
Check that the mappings worked right.
- check the sizes of the bam files
- use bamfilesizes to do this.
- check the sizes of the error files.
- if bam files are noticeably larger or smaller, check on those - was the mapping successful?
If those smaller ones are broken, you'll need to re-run the smaller ones by copying the names of only those that failed into the readFileList_test.txt, and running a version of run_bbmap that calls readFileList_test.txt.

##**Step 9: Merge the read files into 1 bam file per genome**  
Samtools is installed on CHTC. Therefore, you can cd to each metageome directory and run the following, editing the name of the metagenome within each directory. So I used variants of the following:  

```
cd IHXX
ls # to check there is no merged file in there already
time samtools merge all_IHXX.bam *.bam
```

Each merging takes maybe 5 min. Of course, you could do this iteratively, running each mapping directly after one another using a script. However, I would not recommend this because that could tie up the Condor scheduler for a long time, so I suggest breaking up the process and running these lines individually.  



##**Step 10: Move the merged bam files from their individual folders**  

Move the merged bam files from their home folders to a single folder together.  

Run:  

```
./move_merged_bam.sh
```

This moves files to merged_bam.  

Run the following:  

```
./gzip_bam.sh
```

This gzips bam files in merged_bam.  

##**Step 11: Import bamfiles to Zissou and count reads per gene** 

Run:

nohup ./receive_gz.sh &

This code:  
- retrieves a gzipped merged mapping file from CHTC
- stores the file temporarily in Zissou
- unzips the gzipped file
- converts bam to sam
- uses sam file as input into htseq-count, which counts the number of reads mapping to each gene.
- **NOTE!** The input to this script is whatever metagenome mapping files you've listed in readFileList_test.txt. Note this is not the same as the readFileList_test that you used to make run_bbmap.sub work; the difference is that only the 4-letter or 11-letter metagenome code is listed in the text file, not both the code and the names of the metagenome subset files, such as IHPNaa.
- **NOTE!** Listing all metagenome files in readFileList_test causes problems on Zissou. I recommend running ~8 at a time for maximum effiency while (hopefully) avoiding errors.
  
##**Step 12: Quantify reads mapped per COG**

This section is under construction. It will describe how to the script **normalize_counts.R** transforms the count files output by htseq-count into a table of normalized COG abundances. Rows of this table are metagenome sample dates, and columns of this table are normalized COG abundances for that date. **R**eads **P**er **K**ilobase of reference gene per **M**illion mapped reads (RPKM) normalizes the read count data and generates an approximation of the relative number of fragments of a particular gene found in a given sample. Normalizing by RPKM enables comparison of COG abundance both within and across samples. x
  
Additional notes:  

Because the Mendota coassembly is so large (4.2GB), a protocol separate from this one is needed to describe the required testing. There are different. These include:
- Splitting the coassembly into smaller chunks
- Testing for memory when mapping to these smaller chunks.
- Transferring multiple coassembly subsets using scp.
