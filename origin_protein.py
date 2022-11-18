import sys 
import re
def main(pep_file,proteinGroup_file,reffile):
    f1 = open(pep_file,"r")
    f2 = open(proteinGroup_file,"r")
    fo = open(reffile,"a")
    dic_cds = {}
    list1 = []
    for line in f1 :   # cds_file
        line = line.strip()
        if line.startswith('>'):   
            name = line.replace(">","")
        else:
            dic_cds[name] = line
    for line1 in f2:
        genename = line1.split("\t")[6].split(";")
        list1.append(genename)
    for i in list1:
        for k,v in dic_cds.items():
            gene = k.split("_")[0]
            if gene in i:
                fo.write(">"+k+"\n"+v+"\n")
main(sys.argv[1],sys.argv[2],sys.argv[3])
