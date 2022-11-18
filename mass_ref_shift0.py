# Author: shunkai chen 
# Date: 2020-11-24 23:28:24 
# Last Modified by:   shunkai chen 
# Last Modified time: 2020-11-24 23:28:24  
# describe:构建reference
 

import sys 
import re

def main(cds_file,pep_file,shift_file,proteinGroup_file,outfile): 
    f1 = open(cds_file,"r")
    f2 = open(pep_file,"r") 
    f3 = open(shift_file,"r")
    f4 = open(proteinGroup_file,"r")
    fo = open(outfile,"w") 
    dic_cds = {}
    dic_pep = {}
    detected_gene = []
    for line in f1 :   # cds_file
        line = line.strip()
        if line.startswith('>'):    #判断字符串是否以‘>开始’
            name=line.replace(">","").split("_")    #以空格为分隔符，并取序列为0的项。
            name = name[0]+"_"+name[1]+"_"+name[2]
            dic_cds[name]=''
        else:
            dic_cds[name] = line
    #print(dic_cds)

    for lin in f2: # pep_file,use to finder wether have the same sequence in the file
        lin = lin.strip()
        if lin.startswith(">"):
            seq = lin.replace(">","").split("_")
            gene_name = seq[0]+"_"+seq[1]+"_"+seq[2]
            dic_pep[gene_name] = ""
        else:
            dic_pep[gene_name] += lin.replace("\n","")
    
    for l in f4:            #构建已有的基因的序列
        detect_pep = l.strip().split("\t")[6]
        # print(detect_pep)
        detected_gene.append(detect_pep)
    for li in f3:       # 利用文件构建正负移码的蛋白质序列
        li = li.strip().split(",")
        seq1 = li[0].split("_")
        trs_name = seq1[0]+"_"+seq1[1]+"_"+seq1[2]
        shift_motif = li[1].split("_")[1]
        shift_pep_motif = li[1].split("_")[0]
        direct = li[-2]
        gene_we_have = trs_name.split("_")[0]
        for j in range(0,len(shift_pep_motif)):
            shift_pep_site = int(seq1[-1])+int(j)*3
            shift_site = int((shift_pep_site-100)/3)
       
            shift_minus_1_cds = dic_cds[trs_name][shift_pep_site-1:-2]
            shift_plus_1_cds = dic_cds[trs_name][shift_pep_site+1:-3]
        
            down_minus1_pep = translate_dna(shift_minus_1_cds).split("*")[0]   
            down_plus1_pep = translate_dna(shift_plus_1_cds).split("*")[0]

            pepupstream = dic_pep[trs_name][0:shift_site]
            up_by_R ="R" + (pepupstream.split("R")[-1])
            up_by_k ="K" + (pepupstream.split("K")[-1])

            if len(up_by_R) > len(up_by_k) :
                set_sw_minus = "off"
                set_sw_plus = "off"
                minus_pep = up_by_k + down_minus1_pep
                plus_pep = up_by_k + down_plus1_pep
                for K,v in dic_pep.items():
                    if  minus_pep in v :
                        set_sw_minus = "on"
                    if plus_pep in v: 
                        set_sw_plus = "on"
                    else:
                        continue
                if gene_we_have in detected_gene :
                    if direct == "1":
                        if set_sw_plus == "off" and len(plus_pep) >6 and down_plus1_pep != "*" :
                            fo.write(">"+trs_name+":"+"plus"+":"+shift_motif+":"+str(shift_site)+":"+"shift"+str(j)+"\n"+plus_pep+"\n")
                    if direct == "2":                
                        if set_sw_minus == "off" and len(minus_pep) > 6 and down_minus1_pep != "*" :
                            fo.write(">"+trs_name+":"+"minus"+":"+shift_motif+":"+str(shift_site)+":"+"shift"+str(j)+"\n"+minus_pep+"\n")
                else:
                    continue
 
            if len(up_by_R) < len(up_by_k):
                set_sw_minus = "off"
                set_sw_plus = "off"
                minus_pep = up_by_R + down_minus1_pep
                plus_pep = up_by_R + down_plus1_pep   
                for K,v in dic_pep.items():
                    if  minus_pep in v :
                        set_sw_minus = "on"
                    if  plus_pep in v :
                        set_sw_plus = "on"
                    else:
                        continue

                if gene_we_have in detected_gene :
                    if direct == "1":
                        if set_sw_plus == "off" and len(plus_pep) >6 and down_plus1_pep != "*" :
                            fo.write(">"+trs_name+":"+"plus"+":"+shift_motif+":"+str(shift_site)+":"+"shift"+str(j)+"\n"+plus_pep+"\n")
                    if direct == "2":
                        if set_sw_minus == "off" and len(minus_pep) > 6 and down_minus1_pep != "*":
                            fo.write(">"+trs_name+":"+"minus"+":"+shift_motif+":"+str(shift_site)+":"+"shift"+str(j)+"\n"+minus_pep+"\n")
                else:
                    continue
            else:
                continue
        
    fo.close()


def cut(obj, sec):
    return [obj[i:i+sec] for i in range(0,len(obj),sec)]

def translate_dna(sequence):
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }
    codons = [] # Create a codon list to store codons generated from coding seq.
    for i in range(len(sequence)//3):
        if sequence[3*i:3*i+3] in codontable.keys():
            codons.append(sequence[3*i:3*i+3])
    protein_sequence = ''.join([codontable[codon] for codon in codons]) #Translate condons to protein seq.
    return(protein_sequence) 

 
main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])