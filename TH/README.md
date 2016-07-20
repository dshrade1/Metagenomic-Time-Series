## Introduction

This is a protocol to quantify changes in gene content over time based on a time series of metagenomic reads. 

The main steps, generally, include the following: 

1. *Map* a time series of metagenomic reads to a combined assembly of the time series of metagenomes from that sample site.  
2. *Count* the number of reads that align to each locus tag defined in the combined assembly.  
3. *Normalize* the read counts per locus tag per date by reads per kilobase of sequence per million mapped reads (RPKM)
4. *Annotate* the locus tags to quantify the normalized quantity of each COG, KO term, Pfam, etc. over time.

