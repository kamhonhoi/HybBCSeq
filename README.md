## Overview

HybBCSeq is a suite of bioinformatics tools that is used to process and analyze next-generation sequening data generated by the Hybridoma Barcoded Sequencing workflow.

- Demultiplex and bin NGS data to its well of origin
- Cleaning of the NGS data to screen for productive antibody variable domain sequences
- Report the antibody varaible domain sequences for each well

## Prerequisite

- Ubuntu 16.04
- Puthon 2.7
- virtualenv 15.1.0
- Git 2.7.4
- flash 1.2.11 (Included in the HybBCSeq-working directory)

## Installation/Download
- Perform git clone with the following command:
```
git clone https://github.com/kamhonhoi/HybBCSeq.git
```

## Usage

1.  Retrieve raw NGS sequence files (.gz extension) from the source sequencer location to the HybBCSeq-working/samples directory

2.  In order to run the provided scripts, activate virtualenv with the following command:
```
source HybBCSeq-venv/bin/activate

    - Note: to end virtualenv session, use command --- deactivate
    - Change into the HybBCSeq-working directory (i.e. cd HybBCSeq-working)


3.	Merge pair-end reads and re-label output files with desired labeling
##### a.	Program used: flash
##### b.	Usage example: 
###### i.	./flash –r 300 –f 500 –s 50  samples/NGS-R1.fastq.gz samples/NGS-R2.fastq.gz –o samples/NGS-merged
##### c.	Arguments explained:
###### i.	–r : sequence read length per read direction (for MiSeq 2x300, set read length to 300)
###### ii.	–f : expected merged read fragment length
###### iii.	–s : standard deviation from expected read fragment length
###### iv.	Locations of the NGS R1 and R2 sequence files
###### v.	–o : output location and custom prefix
##### d.	Outputs: please refer to the flash help for explanations on the generated files; in particular, the file with .extendedFrags.fastq extension is the merged file needed for next step 

4.	Demultiplexing the merged sequences to wells
a.	Script used: BarcodedSeq-demultiplex.py
b.	Usage example:
i.	python BarcodedSeq-demultiplex.py barcodes.fna samples/NGS-merged.extendedFrags.fastq
c.	Arguments explained:
i.	barcodes.fna : the FASTA file containing the corresponding barcodes for Row and Columns (i.e. VH_barcodes.fna or VK_barcodes.fna)
ii.	merged_sequence file : location of the merged sequence file
d.	Outputs: -demux.csv is the file needed for next step. –demux.log is the log file for the demultiplexing process. –demux-unfoundBC.csv is the file containing sequences without detectable barcodes. –unfoundflag.fna is the FASTA file containing sequences without detectable flag.
5.	Cleaning up the demultiplex sequences
a.	Script used: BarcodedSeq-cleanup.py
b.	Usage example:
i.	python BarcodedSeq-cleanup.py  motif.MotifT  samples/demux.csv
c.	Arguments explained:
i.	motif.MotifT : Location of the probability table flanking the mouse variable domain
ii.	demux.csv : Location of the demultiplexed CSV file
d.	Outputs: -cleaned.csv is the cleaned file for the next step
6.	Consolidating cleaned demultiplexed sequences
a.	Script used: BarcodedSeq-consolidate.py
b.	Usage example:
i.	python BarcodedSeq-consolidate.py –mr 2 –ml 300 samples/cleaned.csv
c.	Arguments explained:
i.	-mr: minimum read counts to be considered for subsequent analysis
ii.	cleaned.csv: Location of the cleaned multiplexed file
d.	Outputs: -cons.csv is the file needed for the next step; -cons-parametersLog.txt contains the arguments parameters used
7.	Reporting the representative sequences for each well
a.	Script used: BarcodedSeq-report.py
b.	Usage example:
i.	python BarcodedSeq-report.py –n 20 samples/cons.csv
c.	Arguments explained:
i.	–n : the number of iterations; higher number increase yields at the expense of representative sequence quality
ii.	cons.csv : Location of the consolidated file
d.	Outputs: -report.csv is the final report file; -report.log reports the number of wells reported
```




## Citation
