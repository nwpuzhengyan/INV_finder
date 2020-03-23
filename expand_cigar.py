import argparse
from PIL import Image
import numpy as np
from pyfaidx import Fasta

#the function to expand cigar
def expand_cigar(cigar,reverse,pos,genes,read_str):
    cigar_str=''
    i=0
    S_clipping=[0,0]
    sign=0
    while i<len(cigar):
        num=''
        while cigar[i].isdigit():
            num=num+cigar[i]
            i+=1
        if sign==0 and cigar[i]=='S':
            S_clipping[0]=int(num)
        if sign==1 and cigar[i]=='S':
            S_clipping[1] = int(num)
        sign=1
        if (cigar[i]!='S' and cigar[i]!='H'):
            cigar_str=cigar_str+cigar[i]*int(num)
        i+=1
    ref=genes[chr][pos:pos+len(cigar_str)]
    read=read_str[S_clipping[0]:len(read_str)-S_clipping[1]]
    cigar_str2=''
    pos=0
    pos1=0
    pos2=0
    while pos<len(cigar_str) and pos1<len(read) and pos2<len(ref):
        if cigar_str[pos]=='M':
            if str(ref[pos2]).upper()==read[pos1].upper():
                cigar_str2+='M'
            else:
                cigar_str2 += 'X'
            pos1 += 1
            pos2 += 1
        elif cigar_str[pos]=='D':
            cigar_str2 += 'D'
            pos2+=1
        else:
            pos1+=1
        pos+=1
    if reverse=='False':
        return cigar_str2
    else:
        return cigar_str2.lower()
#the function to merge cigar of same read
def merge_cigar(cigar1, cigar2, pos1, pos2):
    if pos1<pos2:
        n=pos2-pos1-len(cigar1)
        cigar_str=cigar1+'0'*n+cigar2
        res=[pos1,cigar_str]
    else:
        n = pos1 - pos2 - len(cigar2)
        cigar_str = cigar2 + '0' * n + cigar1
        res = [pos2, cigar_str]
    return res
#the function to read the file
def read_cigrarfile(TXT,fasta,chrsign):
    genes = Fasta(fasta)
    f = open(TXT, "r")
    line = f.readline()
    read_cigar = {}
    while line:
        line = line.split()
        if chrsign!=line[6]:
            line = f.readline()
            continue
        if line[0] in read_cigar:
            pos1 = int(read_cigar[line[0]][0])
            cigar1 = read_cigar[line[0]][1]
            pos2 = int(line[1])
            cigar2 = expand_cigar(line[2], line[4], int(line[1]), genes, line[5])
            read_cigar[line[0]] = merge_cigar(cigar1, cigar2, pos1, pos2)
        else:
            read_cigar[line[0]] = [line[1], expand_cigar(line[2], line[4], int(line[1]), genes, line[5])]
        line = f.readline()
    f.close()
    return read_cigar
#the function to create a matrix
def create_matrix(read_cigar,total_read):
    matrix = np.zeros((len(total_read), matrix_len))
    i = 0
    for key in sorted(total_read):
        if (key in read_cigar.keys())==False:
            matrix[i]==0
            i+=1
            continue
        s = int(read_cigar[key][0])
        e = int(read_cigar[key][0]) + len(read_cigar[key][1])
        s = max(0, s - WIN_S)
        e = min(matrix_len, e - WIN_S)
        for pos in range(0, e - s):
            if WIN_S - int(read_cigar[key][0]) > 0:
                clipping = WIN_S - int(read_cigar[key][0])
            else:
                clipping = 0
            if read_cigar[key][1][clipping + pos] == 'D' or read_cigar[key][1][clipping + pos] == 'd':
                matrix[i][s + pos] = 2
            elif read_cigar[key][1][clipping + pos] == 'X' or read_cigar[key][1][clipping + pos] == 'x':
                matrix[i][s + pos] = 4
            elif read_cigar[key][1][clipping + pos].isupper():
                matrix[i][s + pos] = 1
            elif read_cigar[key][1][clipping + pos] == '0':
                matrix[i][s + pos] = 0
            else:
                matrix[i][s + pos] = 3
        i += 1
    return matrix
#the function to create a image
def create_image(matrix,Image_name,n):
    c = Image.new("RGB", (n*matrix.shape[0], Image_len))
    image_s = int(matrix_len / 2) - int(Image_len / 2)
    image_e = int(matrix_len / 2) + int(Image_len / 2)
    for i in range(0, matrix.shape[0]):
        for j in range(image_s, image_e):
            for k in range(0,n):
                if matrix[i, j] == 0:
                    c.putpixel([i*n+k, j - image_s], (0, 0, 0))
                elif matrix[i, j] == 2:
                    c.putpixel([i*n+k, j - image_s], (255, 0, 0))
                elif matrix[i, j] == 3:
                    c.putpixel([i*n+k, j - image_s], (0, 200, 255))
                elif matrix[i, j] == 4:
                    c.putpixel([i*n+k, j - image_s], (0, 255, 0))
                else:
                    c.putpixel([i*n+k, j - image_s], (255, 255, 0))
    c.save(Image_name)
#the function to expand matrix
def expand_matrix(matrix,n):
    matrix2 = np.zeros((n*matrix.shape[0],matrix.shape[1]))
    for i in range(0,matrix.shape[0]):
        for j in range(0,n):
            matrix2[i*n+j,]=matrix[i,]
    return matrix2
#the function to calculate the type of reads
def calculate_type(matrix1,matrix2):
    hg19_type1,hg19_type2,hg19_type3,hg19_type4 = 0,0,0,0
    rev_hg19_type1, rev_hg19_type2, rev_hg19_type3, rev_hg19_type4 = 0, 0, 0, 0
    extra_read=0
    image_s = int(matrix_len / 2) - int(Image_len / 2)
    image_e = int(matrix_len / 2) + int(Image_len / 2)
    for i in range(0, matrix1.shape[0]):
        plus1 = 0
        plus2 = 0
        rev = 0
        rev_plus1=0
        rev_plus2 = 0
        rev_rev=0
        for j in range(image_s, image_e):
            if rev==0 and matrix1[i, j]==1:
                plus1=1
            elif rev==1 and matrix1[i, j]==1:
                plus2=1
            if matrix1[i, j]==3:
                rev=1

            if rev_rev==0 and matrix2[i, j]==1:
                rev_plus1=1
            elif rev_rev==1 and matrix2[i, j]==1:
                rev_plus2=1
            if matrix2[i, j]==3:
                rev_rev=1

        if (plus1 == 1 and plus2 == 1 and rev == 1) and (rev_plus1==1 and rev_rev==0):
            hg19_type1 += 1
        elif rev == 1 and (plus1 == 1 or plus2 == 1) and (rev_plus1==1 and rev_rev==0):
            hg19_type2 += 1
        elif plus1 == 1 and rev==0 and (rev_plus1==1 and rev_rev==0):
            hg19_type3 += 1
        elif rev == 1 and plus1==0 and (rev_plus1==1 and rev_rev==0):
            hg19_type4 += 1
        elif (rev_plus1 == 1 and rev_plus2 == 1 and rev_rev == 1) and (plus1==1 and rev==0):
            rev_hg19_type1 += 1
        elif rev_rev == 1 and (plus1 == 1 or plus2 == 1) and (plus1==1 and rev==0):
            rev_hg19_type2 += 1
        elif rev_plus1 == 1 and rev_rev==0 and (plus1==1 and rev==0):
            rev_hg19_type3 += 1
        elif rev_rev == 1 and rev_plus1==0 and (plus1==1 and rev==0):
            rev_hg19_type4 += 1
        else:
            extra_read+=1
    support_read=hg19_type1+hg19_type2+hg19_type3
    unsupport_read=rev_hg19_type1+rev_hg19_type2+rev_hg19_type3
    return(str(matrix1.shape[0]) + ' ' + str(hg19_type1) + ' ' + str(hg19_type2) + ' ' + str(hg19_type3) + ' ' + str(hg19_type4)+\
           ' ' + str(rev_hg19_type1) + ' ' + str(rev_hg19_type2) + ' ' + str(rev_hg19_type3) + ' ' + str(rev_hg19_type4)+\
           ' '+str(extra_read)+' '+str(support_read)+' '+str(unsupport_read))

parser = argparse.ArgumentParser()
parser.add_argument("TXT1", help='cigar file')
parser.add_argument("fasta1", help='fasta file')
parser.add_argument("TXT2", help='cigar file')
parser.add_argument("fasta2", help='fasta file')
parser.add_argument("INV", help='fasta file')
args = parser.parse_args()

read_type= open('read_type.txt', 'w')
f = open(args.INV, "r")
line = f.readline()
while line:
    line = line.split()
    read_type.write(str(line[0])+':'+str(line[1])+'-'+str(line[2])+' ')
    chr = line[0]
    pos1 = int(line[1])
    pos2=int(line[2])
    chrsign=chr+':'+str(pos1)
    INV_len=pos2-pos1
    matrix_len = int(INV_len*3)
    Image_len = int(INV_len*3)
    WIN_S = pos1 - int(matrix_len / 2)
    time = int(INV_len/2000+1)
    Image_name1 = chr+str(pos1)+"Sniffles_sample.png"
    Image_name2 = chr+str(pos1)+"Sniffles_reference.png"

    read_cigar1 = read_cigrarfile(args.TXT1, args.fasta1,chrsign)
    read_cigar2 = read_cigrarfile(args.TXT2, args.fasta2,chrsign)
    total_read = list(set(sorted(read_cigar1.keys()) + sorted(read_cigar2.keys())))
    matrix1 = create_matrix(read_cigar1, total_read)
    #create_image(matrix1, Image_name1, time)
    matrix2 = create_matrix(read_cigar2, total_read)
    #create_image(matrix2, Image_name2, time)
    type_num=calculate_type(matrix1,matrix2)
    read_type.write(type_num + '\n')

    line = f.readline()
f.close()
read_type.close()