import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("bam", help='bam file')
parser.add_argument("INV", help='INV file')
args = parser.parse_args()

samfile=pysam.AlignmentFile(args.bam, 'rb')
inv_cigar=open(r'Sniffles_cigar.txt','w')
f = open(args.INV, "r")
line = f.readline()
while line:
    line = line.split()
    chr=line[0]
    pos1=int(line[1])
    pos2=int(line[2])
    for read in samfile.fetch(chr,pos1-200,pos2+200):
        if (str(read.query_sequence) != 'None'):
            inv_cigar.write(str(read.query_name) + ' ' + str(read.reference_start) + ' ' + str(read.cigarstring) + ' ' + str(read.tags[2][1]) + '  ' + \
                str(read.is_reverse) + '  ' + str(read.query_sequence)+ '  ' + str(chr)+':'+str(pos1)+ '\n')
    line = f.readline()
samfile.close()
f.close()