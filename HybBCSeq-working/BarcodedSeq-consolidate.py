#!/usr/bin/python

# Kam Hon Hoi
# 08/07/2017

# Purpose: Consolidate sequences further after the cleanup step
#
# Usage: python NintySix-barcode-consolidate.py -h

from argparse import ArgumentParser
from collections import defaultdict
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Arguments handler
parser=ArgumentParser()
parser.add_argument("infile",type=str,help="Provide the cleaned result from the merged data")
parser.add_argument("-mr","--minread",type=int,default=2,help="Minimum number of reads before being considered in analysis. Default: 2 (i.e. eliminating singleton)")
parser.add_argument("-ml","--minlength",type=int,default=300,help="Minimum length required to be considered for analysis. Default: 300")
args=parser.parse_args()

minread=int(args.minread)
minlength=int(args.minlength)

# Extracting filename
filename=args.infile[:-12]

# Outputting a log recording the parameters used
with open("%s-cons-parametersLog.txt"%filename,"w") as fout:
	print >>fout,"Read abundance cutoff: %d"%args.minread
	print >>fout,"Minimum read length: %d"%args.minlength

# Variables initializations
well=defaultdict(list)
wellabundlist=defaultdict(list)
uniq=defaultdict(int)
aaseqtot=defaultdict(int)
aaseqcontam=defaultdict(list)

####----------------------------------------------------------------------------------####
# Steps for catching redundant reads that were not obvious prior to the cleaning step
# Processing the cleaned merged CSV file and recording the sequences into a list
# Note: this unique collapsing is necessary because after cleaning we might have 
#       identical Variable regions that need to be collapsed
####----------------------------------------------------------------------------------####
with open(args.infile,"rU") as f:
	for line in f:
		line=line.strip()
		if 'Well' in line:
			continue
		tmp=line.split(',')
		aaseq=str(Seq(tmp[6],generic_dna).translate())
		if '*' in aaseq or int(tmp[3])<minread or len(tmp[6])<minlength or 'I' in tmp[0]:
			continue
		uniq[tmp[0]+','+tmp[6]]+=int(tmp[3]) # Unique NT sequence count per well
		aaseqtot[aaseq]+=int(tmp[3])
		aaseqcontam[aaseq].append(tmp[0])
	for i in uniq:
		tmp=i.split(',')
		well[tmp[0]].append(tmp[1]) # Storing all unique NT sequences from each well into this dictionary
		wellabundlist[tmp[0]].append(str(uniq[i])+','+tmp[1]) # Storing per well Read Count for the unique NT sequences

# Processing and generating the consensus output for the report or alignment script
with open(filename+'-cons.csv',"w") as fout:
	print >>fout,"Well,Count,Sequence,AASequence,Dominance,Contamination"
	for w in wellabundlist:
		Seq2Cnt=defaultdict(int)
		for s in wellabundlist[w]:
			tmp=s.split(',')
			Seq2Cnt[tmp[1]]=int(tmp[0])
		for i in sorted(Seq2Cnt,key=Seq2Cnt.get,reverse=True):
			aaseq=str(Seq(i,generic_dna).translate())
			print >>fout,"%s,%d,%s,%s,%.4f,%s"%(w,Seq2Cnt[i],i,aaseq,Seq2Cnt[i]*1.0/aaseqtot[aaseq],':'.join(sorted(set(aaseqcontam[aaseq]))))
