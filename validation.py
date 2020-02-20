import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("INV", help='inversion txt')
parser.add_argument("bam", help='G+G bam file')
args = parser.parse_args()

read=pysam.AlignmentFile(args.bam, 'rb')
f = open(args.INV, "r")
line= f.readline()
read_list=list(read)

val_inv = open(r'validated_INV.txt','w')
support=open(r'support_read.txt','w')
unsupport=open(r'unsupport_read.txt','w')
read_num=open(r'read_num.txt','w')
read_MAPQ=open(r'read_MAPQ.txt','w')

while line:
    print (line)
    line=line.split()
    INV_len=int(line[2])-int(line[1])
    val_inv.write(str(line[0]) + ':' + str(line[1]) + '-' + str(line[2]) + ' ')
    read_num.write(str(line[0])+':'+str(line[1])+'-'+str(line[2])+' ')
    read_score={}
    chr_name = str(line[0]) + '_INV'
    cover1=0
    cover2=0
    cover_read1={}
    cover_read2 = {}
    for r in read.fetch(str(line[0]), int(line[1]), int(line[1]) + 1):
        if cover_read1.has_key(r.query_name) == False:
            cover_read1[r.query_name]=1
            cover1+=1
    for r in read.fetch(chr_name, int(line[1]), int(line[1]) + 1):
        if cover_read1.has_key(r.query_name) == False:
            cover_read1[r.query_name]=1
            cover1+=1
    for r in read.fetch(str(line[0]), int(line[2]), int(line[2]) + 1):
        if cover_read2.has_key(r.query_name) == False:
            cover_read2[r.query_name]=1
            cover2+=1
    for r in read.fetch(chr_name, int(line[1]), int(line[1]) + 1):
        if cover_read2.has_key(r.query_name) == False:
            cover_read2[r.query_name]=1
            cover2+=1
    for r in read_list:
        if str(r.reference_name)==str(line[0]) or str(r.reference_name)==chr_name:
            pos1 = r.reference_start+1
            pos2 = r.reference_start + r.reference_length
            if (pos1 < int(line[1]) and pos2 >= int(line[1])) or (
                    pos1 <= int(line[2]) and pos2 > int(line[2])) == True:

                 if read_score.has_key(r.query_name)==False:
                     read_score.setdefault(r.query_name,[]).append(0)
                     read_score.setdefault(r.query_name, []).append(0)
                     read_score.setdefault(r.query_name, []).append(0)
                     read_score.setdefault(r.query_name, []).append(0)
                     read_score.setdefault(r.query_name, []).append(0)
                     read_score.setdefault(r.query_name, []).append(0)

                 if str(r.reference_name)==str(line[0]):
                     print(str(r.query_name)+' '+str(r.reference_length)+' '+\
                           str(r.query_length)+'\n')
                     print(r.cigartuples[0])
                     if (r.tags[2][1]>read_score[r.query_name][1]):
                          read_score[r.query_name][0]=r.mapping_quality
                          read_score[r.query_name][1] = r.tags[2][1]
                          read_score[r.query_name][4] = r.reference_length
                          read_score[r.query_name][5] = r.query_length
                 elif str(r.reference_name)==chr_name:
                     if (r.tags[2][1] > read_score[r.query_name][3]):
                          read_score[r.query_name][2] = r.mapping_quality
                          read_score[r.query_name][3] = r.tags[2][1]
    #print (read_score)
    num1=0
    num2=0
    read1=0
    read2=0
    extra_read=0
    dif_read=0
    for key in read_score:
        read_MAPQ.write(key+' '+str(read_score[key][0])+' '+str(read_score[key][2])+' '+\
                        str(read_score[key][1])+' '+str(read_score[key][3])+ ' ' + str(INV_len)+'\n')
        if (((read_score[key][1]-read_score[key][3])/INV_len)>0.2 or (read_score[key][1]-read_score[key][3]) >50):
            if (abs(read_score[key][4]-read_score[key][5])>0.05*read_score[key][5]):
                dif_read+=1
            read1 += 1
            num1 += 1
            support.write(key+'\n')
        elif ((read_score[key][3]-read_score[key][1])/INV_len)>0.2 or (read_score[key][3]-read_score[key][1]) >50:
            read2 += 1
            num2+=1
            unsupport.write(key+'\n')
        else:
            extra_read+=1

    read_num.write(str(read1) + ' '+str(read2) + ' '+\
                   str(extra_read)+' '+str(cover1)+' '+str(cover2)+' '+str(INV_len)+'\n')
    print (num1)
    print (num2)
    print(dif_read)
    if num1==0 and num2==0:
        rate = 0
    else:
        rate = float(num1) / float(num1 + num2)
    print (rate)

    if (num1<14 and INV_len<=200)  or (num1<28 and INV_len>200) or dif_read>0.5*num1 or extra_read>15 :
        val_inv.write(str(num1) + ' ' + str(num2) + ' ' + str(rate) + ' ' + str(INV_len) + ' ' + 'false_INV' + '\n')
    else:
        if rate>0.7:
            val_inv.write(str(num1)+' '+str(num2)+' '+str(rate)+' '+str(INV_len)+' '+'homo_INV'+'\n')
        elif rate <=0.7 and rate>0.3:
            val_inv.write(str(num1) + ' ' + str(num2) + ' ' + str(rate)+' '+str(INV_len)+' '+'heter_INV' + '\n')
        else:
            val_inv.write(str(num1) + ' ' + str(num2) + ' ' + str(rate) +' '+str(INV_len)+' '+'false_INV'+ '\n')

    line = f.readline()

f.close()
read.close()