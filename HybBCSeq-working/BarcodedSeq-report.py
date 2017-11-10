#!/usr/bin/python

# Kam Hon Hoi
# 08/07/17
#
# Purpose: This script is for applying the recursively dynamic algorithm to report the sequence for each well
#
# Usage: python NintySix-barcode-report.py -h

from argparse import ArgumentParser
from collections import defaultdict
import numpy as np

# Arguments handler
parser=ArgumentParser()
parser.add_argument("query",type=str,help="The consolidated csv file")
parser.add_argument("-n","--niter",type=int,default=20,help="Number of recursive iterations to encourage more yields per well. Default: 10")
args=parser.parse_args()

# Extract filename
filename=args.query[:-9]

# Read cutoff and dominance cutoff functions
def readCutoffCalc(currentIter,totalIter,initcutoff):
	# Exponential decay of the read cutoffs
	return initcutoff*np.exp(-6.5*currentIter/totalIter)
def domCutoffCalc(currentIter,totalIter):
	# AA Sequence dominance with initial value set at 0.05
	dcalc = 0.05*1.4**(currentIter*1.0/totalIter*10)
	if dcalc > 1.0:
		return 1.0
	else:
		return dcalc

# Processing the files to gather info on contamination, and sequence plate dominance information
passnum=0 # To record how many wells have been reported
storage=defaultdict(list) # To store the output
tempstorage=defaultdict(list) # To store output without the good representative entry
goodwell=[] # To record the well that already reported a representative sequence entry
aaseqcontam=defaultdict(list) # For storing contaminating wells
aaseqtot=defaultdict(int) # For storing the AAseq total
tempaaseqcontam=defaultdict(list)
tempaaseqtot=defaultdict(int)
dupid=defaultdict(int) # Dictionary for storing duplicate grouping ID
dupidcnt=0 # Index for duplicate grouping ID

# Input file processing and seiving for the representative sequence
availwellnum=0
seqcntlist=[]
with open(args.query,"rU") as f, open(filename+'-report.csv',"w") as fout:
	# Going through the input file and storing the information in storage
	for line in f:
		if 'Count' in line:
			print >>fout,"Well,Count,Seq,AASeq,Init-Dominance,Init-Purity,Iteration,Grade,Dominance,Purity,DupID"
			continue
		line=line.strip()
		tmp=line.split(',')
		wellnum=tmp[0]
		seqcntlist.append(int(tmp[1]))
		# Tracking the number of unique wells reported in the input file
		if wellnum not in storage:
			availwellnum+=1
		storage[wellnum].append(line)

	# Calculate an initial read count cutoffs by finding the ratio of (ideal avg reads / actual avg reads) ~ SignalToNoise Ratio  
	seqcntlist=np.array(seqcntlist)
	initreadcutoff=(seqcntlist.sum()*1.0/availwellnum)/seqcntlist.mean()
	print "Initial recommended read cutoff: %.2f"%initreadcutoff

	# Recursively and dynamically identify potential good representative clones
	for n in range(args.niter):
		print "Processing iteration %d ..."%(n+1)
		readcutoff=readCutoffCalc(n,args.niter,initreadcutoff)
		print readcutoff
		dominancecutoff=domCutoffCalc(n,args.niter)
		if n==0:
			print dominancecutoff
			for w in storage:
				for i in storage[w]:
					# Ranking for the grade: AAA would be the most confident
					grade=['X','X','X']
					tmp=i.split(",")
					wellnum=tmp[0]
					cnt=int(tmp[1])
					qseq=tmp[2]
					aaqseq=tmp[3]
					dominance=float(tmp[4])
					contam=tmp[5]
					if cnt >= readcutoff: # Reads more than the initial cutoffs which is when n=0, the SNR
						grade[0]='A'
					if dominance >= dominancecutoff: 
						grade[1]='A'
					if contam.count(':')==0: # AA sequence not shared in any other well
						grade[2]='A'
					gradestring=''.join(grade)
					if gradestring=='AAA' or gradestring=='AAX':
						if w not in goodwell:
							if aaqseq not in dupid:
								dupidcnt+=1
								dupid[aaqseq]=dupidcnt
							print >>fout,"%s,%d,%s,%.4f,%s,%d"%(i,n+1,''.join(grade),dominance,contam,dupid[aaqseq])
							goodwell.append(w)
							passnum+=1
						else:
							tempstorage[w].append(i)
							tempaaseqcontam[aaqseq].append(w)
							tempaaseqtot[aaqseq]+=cnt
					else:
						tempstorage[w].append(i)
						tempaaseqcontam[aaqseq].append(w)
						tempaaseqtot[aaqseq]+=cnt
		
		else:
			print dominancecutoff
			storage=tempstorage
			tempstorage=defaultdict(list)
			aaseqcontam=tempaaseqcontam
			aaseqtot=tempaaseqtot
			tempaaseqcontam=defaultdict(list)
			tempaaseqtot=defaultdict(int)
			for w in storage:
				for i in storage[w]:
					# Ranking for the grade: AAA would be the most confident
					grade=['X','X','X']
					tmp=i.split(",")
					wellnum=tmp[0]
					cnt=int(tmp[1])
					qseq=tmp[2]
					aaqseq=tmp[3]
					dominance=cnt*1.0/aaseqtot[aaqseq]
					contam=':'.join(sorted(set(aaseqcontam[aaqseq])))
					if cnt >= readcutoff: # Recursively lower the read count prorated by the number of iterations until it reaches 0
						grade[0]='A'
					if dominance >= dominancecutoff: 
						grade[1]='A'
					if contam.count(':')==0: # AA sequence not shared in any other well
						grade[2]='A'
					gradestring=''.join(grade)
					if gradestring=='AAA' or gradestring=='AAX':
						if w not in goodwell:
							if aaqseq not in dupid:
								dupidcnt+=1
								dupid[aaqseq]=dupidcnt
							print >>fout,"%s,%d,%s,%.4f,%s,%d"%(i,n+1,''.join(grade),dominance,contam,dupid[aaqseq])
							goodwell.append(w)
							passnum+=1
						else:
							tempstorage[w].append(i)
							tempaaseqcontam[aaqseq].append(w)
							tempaaseqtot[aaqseq]+=cnt
					else:
						tempstorage[w].append(i)
						tempaaseqcontam[aaqseq].append(w)
						tempaaseqtot[aaqseq]+=cnt

# Reporting the numbers of well that has a representative sequence
with open(filename+'-report.log',"w") as fout:
	print >>fout,"Number of NGS recovered well: %d"%passnum
