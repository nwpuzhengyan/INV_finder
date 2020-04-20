import pysam
import argparse
def check_cigar(start_pos,cigar,pos1,pos2):
    cigar_str = ''
    i = 0
    S_clipping = [0, 0]
    sign = 0
    while i < len(cigar):
        num = ''
        while cigar[i].isdigit():
            num = num + cigar[i]
            i += 1
        if sign == 0 and cigar[i] == 'S':
            S_clipping[0] = int(num)
        if sign == 1 and cigar[i] == 'S':
            S_clipping[1] = int(num)
        sign = 1
        if (cigar[i] != 'S' and cigar[i] != 'H'and cigar[i]!= 'I'):
            cigar_str = cigar_str + cigar[i] * int(num)
        i += 1
    cigar_pos1=max(pos1-start_pos,0)
    cigar_pos2=min(len(cigar_str),pos2-start_pos)
    inv_region=cigar_str[cigar_pos1-30:cigar_pos2+30]
    Del_len=0
    for i in range(0,len(inv_region)):
        if (inv_region[i]=='D'):
            Del_len+=1
    if float(Del_len)/float(pos2-pos1)>0.5:
        return False
    else:
        return True

def read_bam(samfile,chr,inv_start_pos,inv_end_pos):
    read_list = {}
    for read in samfile.fetch(chr, inv_start_pos - 200, inv_end_pos + 200):
        if read.query_name in read_list:
            ref_start_pos = min(read.reference_start, read_list[read.query_name][0])
            ref_end_pos = max(read.reference_end, read_list[read.query_name][1])
            align_score = read.tags[2][1] + read_list[read.query_name][2]
            if read.reference_start > read_list[read.query_name][0]:
                align_dir = read_list[read.query_name][3]+[str(read.is_reverse)]
            else:
                align_dir = [str(read.is_reverse)] + read_list[read.query_name][3]
            if read.reference_start < read_list[read.query_name][0]:
                n = read_list[read.query_name][0] - read.reference_end
                read_cigar = read.cigarstring + str(n)+'0'+read_list[read.query_name][4]
            else:
                n = read.reference_start-read_list[read.query_name][1]
                read_cigar = read_list[read.query_name][4]+ str(n)+'0'+read.cigarstring
            read_list[read.query_name] = [ref_start_pos, ref_end_pos, align_score, align_dir,read_cigar]
        else:
            read_list[read.query_name] = [read.reference_start, read.reference_end, read.tags[2][1],
                                          [str(read.is_reverse)],read.cigarstring]
    return read_list

def calculate_read_type(total_read,read_list,rev_read_list,inv_start_pos,inv_end_pos,INV_len):
    hg19_type1, hg19_type2, hg19_type3, hg19_type4, hg19_type5,hg19_type6= 0, 0, 0, 0,0,0
    rev_hg19_type1, rev_hg19_type2, rev_hg19_type3, rev_hg19_type4,rev_hg19_type5,rev_hg19_type6 = 0, 0, 0, 0,0,0
    extra_read=0
    for read_name in sorted(total_read):
        read_overlap = 0
        rev_read_overlap = 0
        if (read_name in read_list.keys())==True:
            if (inv_start_pos > read_list[read_name][1] or inv_end_pos < read_list[read_name][0])==False:
                 read_overlap = min(inv_end_pos, read_list[read_name][1]) - max(inv_start_pos, read_list[read_name][0])
                 rest_len=read_list[read_name][1]-read_list[read_name][0]-read_overlap
        if (read_name in rev_read_list.keys()) == True:
            if (inv_start_pos > rev_read_list[read_name][1] or inv_end_pos < rev_read_list[read_name][0])==False:
                 rev_read_overlap = min(inv_end_pos, rev_read_list[read_name][1]) - max(inv_start_pos,rev_read_list[read_name][0])
                 rev_rest_len = rev_read_list[read_name][1] - rev_read_list[read_name][0] - rev_read_overlap
        if (read_name in read_list.keys()) == False:
            if rev_read_overlap>=10 and rev_rest_len>=10:
                hg19_type5+=1
            else:
                extra_read+=1
            continue
        elif (read_name in rev_read_list.keys()) == False:
            if read_overlap>=10 and rest_len>=10:
                rev_hg19_type5 += 1
            else:
                extra_read+=1
            continue
        else:
            if read_list[read_name][3]==['False','True','False'] and rev_read_list[read_name][3]==['False'] and\
                    rev_read_overlap>=10 and rev_rest_len>=10:
                hg19_type1+=1
            elif rev_read_list[read_name][3]==['False','True','False'] and read_list[read_name][3]==['False'] and\
                    read_overlap>=10 and rest_len>=10:
                rev_hg19_type1+=1
            elif (read_list[read_name][3]==['False','True'] or read_list[read_name][3]==['True','False']) and\
                    rev_read_list[read_name][3]==['False'] and rev_read_overlap>=10 and rev_rest_len>=10:
                hg19_type2 += 1
            elif (rev_read_list[read_name][3]==['False','True'] or rev_read_list[read_name][3]==['True','False']) and\
                    read_list[read_name][3]==['False'] and read_overlap>=10 and rest_len>=10:
                rev_hg19_type2 += 1
            elif read_list[read_name][3]==['True'] and rev_read_list[read_name][3]==['False'] and\
                    rev_read_overlap>=10 and rev_rest_len>=10:
                hg19_type4+=1
            elif rev_read_list[read_name][3]==['True'] and read_list[read_name][3]==['False'] and\
                    read_overlap>=10 and rest_len>=10:
                rev_hg19_type4+=1
            elif rev_read_list[read_name][3]==read_list[read_name][3] and len(rev_read_list[read_name][3])==1:
                ref_len = read_list[read_name][1] - read_list[read_name][0]
                rev_ref_len = rev_read_list[read_name][1] - rev_read_list[read_name][0]
                if (rev_ref_len-ref_len)>=20 and rev_read_overlap>=10 and rev_rest_len>=10:
                    if check_cigar(rev_read_list[read_name][0],rev_read_list[read_name][4],inv_start_pos,inv_end_pos):
                        hg19_type3+=1
                    else:
                        extra_read+=1
                elif (ref_len-rev_ref_len)>=20 and read_overlap>=10 and rest_len>=10:
                    if check_cigar(read_list[read_name][0], read_list[read_name][4], inv_start_pos,inv_end_pos):
                        rev_hg19_type3+=1
                    else:
                        extra_read+=1
                else:
                     if (rev_read_list[read_name][2]-read_list[read_name][2])>(0.5*INV_len) and INV_len<=1000 and\
                             rev_read_overlap>=10 and rev_rest_len>=10:
                         if check_cigar(rev_read_list[read_name][0], rev_read_list[read_name][4], inv_start_pos,
                                        inv_end_pos):
                             hg19_type6 += 1
                         else:
                             extra_read += 1
                     elif (read_list[read_name][2]-rev_read_list[read_name][2])> (0.5*INV_len) and INV_len<=1000 and\
                             read_overlap>=10 and rest_len>=10:
                         if check_cigar(read_list[read_name][0], read_list[read_name][4], inv_start_pos, inv_end_pos):
                             rev_hg19_type6 += 1
                         else:
                             extra_read += 1
                     else:
                         extra_read+=1
            else:
                extra_read+=1
    support_read = hg19_type1 + hg19_type2 + hg19_type3+ hg19_type5 + hg19_type6
    unsupport_read = rev_hg19_type1 + rev_hg19_type2 + rev_hg19_type3+ rev_hg19_type5 + rev_hg19_type6
    if (support_read+unsupport_read)==0:
        rate=0
    else:
        rate=float(support_read)/float(support_read+unsupport_read)
    if rate > 0.7:
        INV_type = "Homo_INV"
    elif rate <= 0.7 and rate >= 0.3:
        INV_type = "Heter_INV"
    else:
        INV_type = "False_INV"
    return (str(len(total_read)) + ' ' + str(hg19_type1) + ' ' + str(hg19_type2) + ' ' + str(hg19_type3) + ' ' +\
            str(hg19_type4) +' '  + str(hg19_type5) + ' ' +str(hg19_type6) +' '+\
            str(rev_hg19_type1) + ' ' + str(rev_hg19_type2) + ' ' + str(rev_hg19_type3) + ' ' +\
            str(rev_hg19_type4) +' '+ str(rev_hg19_type5) + ' ' + str(rev_hg19_type6) + ' ' +\
            str(extra_read) + ' ' + str(support_read) + ' ' + str(unsupport_read)+' '+str(rate)+' '+str(INV_type))

parser = argparse.ArgumentParser()
parser.add_argument("bam1", help='hg19 bam file')
parser.add_argument("bam2", help='reversed bam file')
parser.add_argument("INV", help='INV file')
args = parser.parse_args()

samfile=pysam.AlignmentFile(args.bam1, 'rb')
rev_samfile=pysam.AlignmentFile(args.bam2, 'rb')
read_type= open('read_type.txt', 'w')
f = open(args.INV, "r")
line = f.readline()
while line:
    line = line.split()
    print(line)
    read_type.write(str(line[0]) + ':' + str(line[1]) + '-' + str(line[2]) + ' ')
    chr = line[0]
    inv_start_pos = int(line[1])
    inv_end_pos = int(line[2])
    INV_len=inv_end_pos-inv_start_pos
    read_list = read_bam(samfile,chr,inv_start_pos,inv_end_pos)
    rev_read_list=read_bam(rev_samfile,chr,inv_start_pos,inv_end_pos)
    total_read = list(set(sorted(read_list.keys()) + sorted(rev_read_list.keys())))
    type_num = calculate_read_type(total_read,read_list,rev_read_list,inv_start_pos,inv_end_pos,INV_len)
    read_type.write(type_num + '\n')
    line = f.readline()
f.close()
read_type.close()
