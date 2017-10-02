__author__ = 'EN'


###########
## USAGE ##
###########

import sys
import re
import os
import pysam
import getopt



def read_prediction(RNAFile, rna_reads, RiboFile, ribo_reads, inputFile, outputFile):

	samRNA = pysam.Samfile(RNAFile, "rb")
	samRibo = pysam.Samfile(RiboFile, "rb")

	header = ["ORF_locus","strand","length","ORF_type","start_codon","Reference Annotation","Nterminal Support","Ribo_RPKM", "RNA_RPKM", "TE"]
	
	outFile = open(outputFile, 'w')
	outFile.write('\t'.join(header) + '\n')

	inFile = open(inputFile, 'r')
	line = inFile.readline()

	while line != '':
		if not (line.startswith('ORF_locus')):
			fields = line.split()

			orf_id = str(fields[0]) 
			strand = str(fields[1])
			#length = fields[2]
			ORF_type = fields[3]
			start_codon = fields[4]
			ref = fields[5]
			nterm = fields[6]
			#dist = fields[7]

			region = orf_id.split(":")[0]
			position = orf_id.split(":")[1]
			start = int(position.split("-")[0])
			stop = int(position.split("-")[1])

			length = stop - start + 1
			# RNA RPKM
			rna_count = 0
			rna_rpkm = 0.0
			if strand == '+':	# count read on nagative strand
				for read in samRNA.fetch(str(region),int(start),int(stop)):
					if not read.is_reverse:
						rna_count = rna_count + 1

			elif strand == '-':
				for read in samRNA.fetch(str(region),int(start),int(stop)):
					if read.is_reverse:
						rna_count = rna_count + 1

			rna_rpkm = float(rna_count*1000000000)/(length*rna_reads)

			# Ribo RPKM
			ribo_count = 0
			ribo_rpkm = 0.0
			if strand == '+':	# count read on nagative strand
				for read in samRibo.fetch(str(region),int(start),int(stop)):
					if not read.is_reverse:
						ribo_count = ribo_count + 1

			elif strand == '-':
				for read in samRibo.fetch(str(region),int(start),int(stop)):
					if read.is_reverse:
						ribo_count = ribo_count + 1

			ribo_rpkm = float(ribo_count*1000000000)/(length*ribo_reads)

			TE = 0.0
			if rna_rpkm == 0:
				TE = 0.0
			else:
				TE = ribo_rpkm/rna_rpkm

			new_line = fields
			new_line.append(str(ribo_rpkm))
			new_line.append(str(rna_rpkm))
			new_line.append(str(TE))

			outFile.write('\t'.join(new_line) + '\n')
		line = inFile.readline()

	samRNA.close()
	samRibo.close()



if __name__=='__main__':

	RNAFile = sys.argv[1]
	rna_reads = int(sys.argv[2])
	RiboFile = sys.argv[3]
	ribo_reads = int(sys.argv[4])
	inputFile = sys.argv[5]
	outputFile = sys.argv[6]
	read_prediction(RNAFile, rna_reads, RiboFile, ribo_reads, inputFile, outputFile)



