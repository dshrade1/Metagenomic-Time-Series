#!/usr/bin/env python 
import sys
from Bio import SeqIO

inFile  = sys.argv[1]
SeqIO.convert(inFile, "fastq", inFile + ".fasta", "fasta")

#short_sequences = [] # Setup an empty list
#for record in SeqIO.parse(open("cor6_6.fastq", "rU"), "fastq") :
#    if len(record.seq) > 150 :
#        # Add this record to our list
#        short_sequences.append(record)
 
# print "Found %i short sequences" % len(short_sequences)

# http://biopython.org/wiki/SeqIO#File_Format_Conversion
