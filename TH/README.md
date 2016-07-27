## Introduction

This is a protocol to quantify changes in gene content over time based on a time series of metagenomic reads. 

The main steps, generally, include the following: 

1. **Map** a time series of metagenomic reads to a combined assembly of the time series of metagenomes from that sample site.  
2. **Count** the number of reads that align to each locus tag defined in the combined assembly.  
3. **Normalize** the read counts per locus tag per date by reads per kilobase of sequence per million mapped reads (RPKM)
4. **Annotate** the locus tags to quantify the normalized quantity of each COG, KO term, Pfam, etc. over time.

## Requirements

To run this protocol, you'll need a few things:  

* 1\. An account with the Center for High-Throughput Computing for 
    *  High-Throughput Computing node (to perform mapping) with Python (I used 2.7), biopython (I used v1.67), and samtools 1.3 installed. 
    *  Gluster account (to store metagenome data in CHTC)
    *  SQUID account (to store combined assembly in CHTC)
* 2\. A combined assembly of the time series for your sample site of interest 
    *  in .fna file format  
    *  associated annotation files (from JGI, contact Sarah Stevens sstevens2@wisc.edu) 
* 3\. A directory of the time series of metagenomic read files for your sample site of interest 
    *  unassembled, merged reads in fasta or fastq format
    *  additionally, a text file or csv of sampling dates for each metagenome
* 3\. A GFF file of the combined assembly.

## Protocol

The following are the steps of the protocol in detail. Different steps require different computing needs, so the protocol separated by the location where the computing is happening.

### in CHTC
* 1\. Place the required input data in the appropriate locations. Cyberduck is a handy SFTP client for doing this. This means: 
    *  coassembly file (zipped) in your SQUID directory
    *  directory of merged metagenomic reads in your Gluster directory
    *  have the python, samtools, htseq-count, and other software installed in the submit node.
* 2\. Clone the repo Metagenomic-Time-Series-CHTC into your home directory on the CHTC submit node. This contains: 
    *  one directory for each sample site (TH = Trout Bog Hypolimnion, TE = TB Epi, ME = Mendota)
    *  folders containing scripts for preprocessing and mapping specific to that sample site
    *  map.dag script in each directory
* 3\. Update your $PATH in the submit node so Condor looks for the version of Python YOU installed, rather than its default version, first.
```
cd ~/Metagenomic-Time-Series-CHTC/global/
export PATH=$(pwd)/python/bin:$PATH
echo $PATH
```
* 4\. Move to the directory for your sample site of interest.
* 5\. If you wish to adjust mapping parameters, open the .dag file. Edit the value in the "VAR" line(s). For example, to increase the minimum percent identity threshold for mapping with BBMap, change *VAR mapping minid="76"* to *VAR mapping minid="92"*. Also edit the name of your squid directory and the name of the directory in Gluster containing your metagenome reads. *Note to self: make these options available in dag file.*
* 6\. Submit the dag file
``` 
condor_submit_dag map.dag
```
This will execute the following steps:  
    *  **Preprocess metagenomic read files for mapping.** For fasta files, line-wrapping is removed, and each metagenomic read file is split into many <20MB pieces so that they can be transported to CHTC's compute nodes for faster mapping down the line. Additional formatting steps are performed for fastq files.  
    * **Map metagenomic read files to combined assembly.** Once the preprocessing has successfully completed, metagenomic read files are mapped to the coassembly using BBMap software. One directory is created for each metagenome, and .bam output files are sent to each metagenome's mapping directory.
    * **Count the number of reads that mapped to each locus tag in the combined assembly.** The Python module HTSeq-count is called to count the number of reads that mapped to each locus tag whether those are coding sequences, repeat regions, tRNA, or rRNA.

* 2\. Transfer the output directory of reads counts to your local machine. 
  
### on your local machine

* 1\. Clone this GitHub repository to your local machine.  
* 2\.  Run count_features.R. This will read the htseq-count read count files into R and separate out mapped reads from unmapped reads.  
* 3\. Run counts_analysis.R. This will answer some relevant questions about the counts output, generating and saving descriptive plots, such as:
    * What is the overall proportion of reads that mapped as CDS, repeat_region, rRNA, tRNA, etc.?
    * What is the total % of reads mapped in each metagenome?
    * Which dates have the highest % of total reads mapping, and which have the lowest?
    * How many reads mapped as each feature type on each date? 
    * Are the ratios of [CDS + aligned noncoding]:total reads consistent?
    * For samples taken on the same date, how similar are the percents of reads mapped?
* 4\. Run normalize_counts.R This will normalize read counts by RPKM.  
* 5\.  Run normalize_analysis.R This will answer questions such as:  
    * What are the sums of the RPKM values per metagenome?  
    * Are they consistent?
* 6\.  Run compare_annotations.R. This will help you identify which annotation scheme (COG, KO, EC, Pfam, etc.) is most appropriate for the data set you are using.
* 7\.  Run annotate_genes.R. This will allow you to choose an annotation scheme, identify the "best" annotation for each locus tag (first by E-value, then by bit score, then by %identity*alignment_length), then annotate each locus tag with that best annotation. The RPKM values for all genes with a given annotation are summed, so that each, say, COG, has 1 RPKM value per date. This script creates a time series of RPKM values for each, say, COG, in the annotation file.
