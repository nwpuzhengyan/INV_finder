import argparse
from PIL import Image
import numpy as np

def expand_cigar(cigar):
    cigar_str=''
    i=0
    while i<len(cigar):
        num=''
        while cigar[i].isdigit():
            num=num+cigar[i]
            i+=1
        if (cigar[i]!='S' and cigar[i]!='H' and cigar[i]!='I'):
            cigar_str=cigar_str+cigar[i]*int(num)
        i+=1
    return cigar_str


parser = argparse.ArgumentParser()
parser.add_argument("TXT", help='cigar file')
args = parser.parse_args()

f = open(args.TXT, "r")
line= f.readline()
read_cigar={}

while line:
    line = line.split()
    if line[0] in read_cigar:
        if (int(line[3])>int(read_cigar[line[0]][2])):
            read_cigar[line[0]]=[line[1], expand_cigar(line[2]),line[3]]
    else:
        read_cigar[line[0]] = [line[1], expand_cigar(line[2]), line[3]]


    line = f.readline()

f.close()

matrix=np.zeros((len(read_cigar),20000))

i=0
for key in read_cigar:
    s=int(read_cigar[key][0])
    e=int(read_cigar[key][0])+len(read_cigar[key][1])
    WIN_S=100038214
    s=max(0,s-WIN_S)
    e=min(20000,e-WIN_S)
    for pos in range(0,e-s):
        if WIN_S-int(read_cigar[key][0])>0:
            clipping=WIN_S-int(read_cigar[key][0])
        else:
            clipping=0
        if read_cigar[key][1][clipping+pos]=='D':
            matrix[i][s+pos] = 2
        else:
            matrix[i][s + pos] = 1
    i+=1

c = Image.new("RGB", (len(read_cigar),2000))

image_s=9000
image_e=11000

for i in range (0,len(read_cigar)):
    for j in range (image_s,image_e):
       if matrix[i,j]==0:
           c.putpixel([i,j-image_s],(0,0,0))
       elif matrix[i,j]==2:
           c.putpixel([i, j - image_s], (255, 0, 0))
       else:
           c.putpixel([i, j-image_s], (255, 255, 0))

c.show()
c.save("large_sample.png")

