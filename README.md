## Introduction

This is a protocol to map metagenomic reads to a metagenomic coassembly. Steps covered are as follows:  

1. Set up directory for pre-mapping formatting in lab server (Zissou)  
2. Set up directory for mapping in compute cluster (CHTC)  
3. Transfer coassembly file to CHTC  
4. Format and transfer metagenomic read files to CHTC  
5. Prep CHTC submit node for mapping  
6. Test mapping in CHTC  
7. Run mapping in CHTC  

##**Step 1: In Zissou, set up directory for pre-mapping formatting steps**##

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
- Python
- Perl
- sshpass

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

##**3. Prep CHTC SQUID folder for mapping.** Add the .fna file for the combined assembly (hereafter, "coassembly") from Zissou to your squid folder on CHTC. 

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

##**Step 4. Format files in Zissou and transfer them to CHTC.**

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

**6a. Copy the first two lines of readFileList.txt into a document called readFileList_test.txt.** It should look like this:  
```
H0586_93745 H0586_93745aa  
H0586_93745 H0586_93745ab  
```

Note: Mapping one 16.9-MB read file to 563-MB coassembly takes  
- 4.03 GB memory
- 37 MB disk space output
- 48.5 seconds to run.
So you'll need to request this much memory.


Because the Mendota coassembly is so large (4.2GB), a protocol separate from this one is needed to describe the required testing. There are different. These include:
Splitting the coassembly into smaller chunks
Testing for memory when mapping to these smaller chunks.
Transferring multiple coassembly subsets using scp.

Future:
Include how long things took.
Include memory usage.
