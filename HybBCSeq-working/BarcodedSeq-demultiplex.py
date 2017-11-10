#!/usr/bin/python

# Kam Hon Hoi
# 08/07/2017

# Purpose: To demultiplex FLASH-merged FASTQ files for the Barcoded sequencing workflow
#
# Usage: python NintySix-barcode-demultiplex.py -h

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from collections import defaultdict
import re
import numpy as np

# Arguments Handler
parser=ArgumentParser()
parser.add_argument("barcodes",type=str,help="Location of corresponding barcodes FASTA file")
parser.add_argument("inputfile",type=str,help="Location of merged FASTQ file")
args=parser.parse_args()

# Extracting input filename prefix
filename=args.inputfile[:-6]

# Generating barcode hash table from the barcodes FASTA file
# Note: key is the barcode and value is the corresponding well position
#       generate 1 bp mismatch variants for the well's hash table
forbc=defaultdict(str)
revbc=defaultdict(str)
forflag=''
revflag=''
nts='ATGC'
with open(args.barcodes,"rU") as f:
	for record in SeqIO.parse(f,"fasta"):
		header=record.id
		seq=str(record.seq)
		if 'reverse' in header:
			revflag=seq[-4:]
			revbc[seq[:-4]]=header.split("_")[0]
			for p in range(10):
				bclist=list(seq[:-4])
				allowmis=list(nts.replace(bclist[p],''))
				for j in allowmis:
					bclist[p]=j
					newbc=''.join(bclist)
					revbc[newbc]=header.split("_")[0]
		else:
			forflag=seq[:4]
			forbc[seq[4:]]=header.split("_")[0]
			for p in range(10):
				bclist=list(seq[4:])
				allowmis=list(nts.replace(bclist[p],''))
				for j in allowmis:
					bclist[p]=j
					newbc=''.join(bclist)
					forbc[newbc]=header.split("_")[0]


# Initialize QC variables
count=defaultdict(int) # Note: the key to this dict is well#_ForBC_RevBC_payloadSeq
qualcount=defaultdict(list) # Note: this stores the min phred score and the average phred score
noBCcount=defaultdict(int) # Note: the no BC dict follows similar setup as the count above except only partial or no well# is reported
noBCqualcount=defaultdict(list)
totseqcount=0
unfoundflag=0
unfoundbc=0
found=0

# Initialize re compile
regex=re.compile(r"%s(.+)%s"%(forflag,revflag))

# Processing the fastq file
with open(args.inputfile,"rU") as f, open(filename+'-unfoundflag.fna',"w") as foutnoflag:
	for record in SeqIO.parse(f,"fastq"):
		totseqcount+=1
		seq=str(record.seq)
		rseq=str(Seq(seq,generic_dna).reverse_complement())
		qual=record.letter_annotations['phred_quality']
		rqual=qual[::-1]
		m=regex.search(seq)
		revm=regex.search(rseq)
		if m:
			seq=m.group(1)
			seqforbc=seq[:10]
			seqrevbc=seq[-10:]
			seq=seq[10:-10]
			qual=qual[m.start()+10:len(seq)]
			k=forbc.get(seqforbc,'Row')+revbc.get(seqrevbc,'Col')+'_'+seqforbc+'_'+seqrevbc+'_'+seq
			if not qual:
				unfoundbc+=1
				noBCcount[k]+=1
				noBCqualcount[k]=[0,0]
				continue
			if seqforbc not in forbc or seqrevbc not in revbc:
				unfoundbc+=1
				noBCcount[k]+=1
				noBCqualcount[k]=[np.min(qual),np.average(qual)]
				continue
			count[k]+=1
			qualcount[k]=[np.min(qual),np.average(qual)]
			found+=1
			continue
		elif revm:
			rseq=revm.group(1)
			seqforbc=rseq[:10]
			seqrevbc=rseq[-10:]
			seq=rseq[10:-10]
			qual=rqual[revm.start()+10:len(seq)]
			k=forbc.get(seqforbc,'Row')+revbc.get(seqrevbc,'Col')+'_'+seqforbc+'_'+seqrevbc+'_'+seq
			if not qual:
				unfoundbc+=1
				noBCcount[k]+=1
				noBCqualcount[k]=[0,0]
				continue
			if seqforbc not in forbc or seqrevbc not in revbc:
				unfoundbc+=1
				noBCcount[k]+=1
				noBCqualcount[k]=[np.min(qual),np.average(qual)]
				continue
			count[k]+=1
			qualcount[k]=[np.min(qual),np.average(qual)]
			found+=1
			continue
		else:
			unfoundflag+=1
			print >>foutnoflag,">%s\n%s"%(record.id,seq)


# Output results and log file
with open(filename+'-demux.csv',"w") as fout:
	print >>fout,"Well,ForBC,RevBC,SeqCount,MinPhred,AvgPhred,Sequence"
	for i in count:
		tmp=i.split('_')
		print >>fout,"%s,%s,%s,%d,%.1f,%.2f,%s"%(tmp[0],tmp[1],tmp[2],count[i],qualcount[i][0],qualcount[i][1],tmp[3])

with open(filename+'-demux-unfoundBC.csv',"w") as fout:
	print >>fout,"Well,ForBC,RevBC,SeqCount,MinPhred,AvgPhred,Sequence"
	for i in noBCcount:
		tmp=i.split('_')
		print >>fout,"%s,%s,%s,%d,%.1f,%.2f,%s"%(tmp[0],tmp[1],tmp[2],noBCcount[i],noBCqualcount[i][0],noBCqualcount[i][1],tmp[3])

with open(filename+'-demux.log',"w") as fout:
	print >>fout,"Total Seqs processed: %d\nSequences with identifiable BCs: %d\nSequences without identifiable BCs: %d\nSequences without flag: %d"%(totseqcount,found,unfoundbc,unfoundflag)
