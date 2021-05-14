#!/usr/bin/python

outFasta = open('Kpneumoniae_mgh78578.cds.fna', 'w')
with open('Klebsiella_pneumoniae_subsp_pneumoniae_mgh_78578.ASM1630v1.cds.all.renamed.format.fa', 'r') as inFasta:
	seq = ''
	firstEntry = inFasta.readline()
	outFasta.write(firstEntry)
	for line in inFasta:
		if line[0] == '>':
			outFasta.write(seq + '\n\n')
			seq = ''
			outFasta.write(line)
		else:
			seq += line.strip()

outFasta.write(seq + '\n')
outFasta.close()
