#!/bin/env python3
import argparse
import os
import sys
import pandas as pd
import numpy as np
from math import ceil
#import logomaker as lm


def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='?',  help='The clusters of peptides ')

    my_parser.add_argument('-db',  nargs='?', default="uniprot_sprot",
                                    help='Blast database ')
    my_parser.add_argument('-evalue', nargs='?', default="10",help='E-value cut off')

    my_parser.add_argument('-op', nargs='?',  default="32767", help='Gap opening penalty')
    my_parser.add_argument('-ep', nargs='?',  default="32767", help='Gap extension penalty')

    return my_parser.parse_args()


def runBlastp (dfTmp, db, evalue, op, ep):
    
    arrTmps = list(dfTmp.index)
    outTmp = open("myTmpPeps.fasta", "w")
    for tmp in arrTmps:
        outTmp.write(">" + tmp + "\n")
        outTmp.write(tmp + "\n")
    outTmp.close()
    
    blastdb = "blastp -db /opt/nextera/galaxy/data/blastdb/" + db
    myCmd1 = blastdb + "  -query myTmpPeps.fasta -num_alignments 1000000000 -num_threads 10 -gapopen " + op + " -gapextend " + ep +" -out myBlastp.out.orig.txt -evalue " + evalue
    os.system(myCmd1)
    myCmd2 = ''' gawk 'length($0) != 0' myBlastp.out.orig.txt > myTmpPeps.out2 '''
    os.system(myCmd2)
    
    ##organize the table
    file0 = open("myTmpPeps.out2")
    arrLine = []
    while True:
        line = file0.readline().rstrip()
        if (len(line) == 0):
            break
        arrLine.append(line)
    file0.close()

    out = open("myTmpPeps.out.txt","w")
    out.write("Peptides" + "\t" + "Accession" + "\t" + "Score(bits)" + "\t" + "E-value" + "\t" + "Identity" + "\t" + "Similarity" + "\t" + "Description" + "\n" )

    for index_i in range(0,len(arrLine)):
        ### "Query=" block
        if(arrLine[index_i].startswith("Query= ")):
            query = arrLine[index_i].split("Query= ")[1].lstrip()
            
            for index_j in range(index_i + 1, len(arrLine)):
                if (arrLine[index_j].startswith("Query= ") ):
                    break
                else:

                    if ("No hits found" in arrLine[index_j ]):
                        out.write(query + "\t" + "No hit found" + "\n")

                    else:

                        ### ">" block
                        if (arrLine[index_j].startswith(">")):
                            pro = ""
                            index_length = 0
                            for index_p in range(index_j, len(arrLine)):
                                if (arrLine[index_p].startswith("Length=")):
                                    index_length = index_p
                                    break
                                else:
                                    pro = pro + arrLine[index_p]

                            pro = pro.replace("> ","")
                            pro_id = pro.split(" ")[0]
                            pro_name = " ".join(pro.split(" ")[1:])

                            score = arrLine[index_length + 1].lstrip().split("Score = ")[1].split(" bits")[0]
                            evalue = arrLine[index_length + 1].lstrip().split("Expect = ")[1].split(", Method")[0]
                            identities = arrLine[index_length + 2].lstrip().split("Identities = ")[1].split(",")[0] 
                            positives = arrLine[index_length + 2].lstrip().split("Positives = ")[1].split(",")[0]

                            out.write(query + "\t" + pro_id + "\t" + score + "\t" + evalue + "\t" + identities + "\t" + positives + "\t" + pro_name + "\n")
            
    out.close()
    
    df_homo = pd.read_csv("myTmpPeps.out.txt",header=0,index_col=0,sep="\t")
    
    #cond1 = df_homo[option + "(%)"] >= threshold
    #cond2 = df_homo["#OverlapAA"] >= num_overlap

    #df_homo = df_homo[cond1 & cond2 ]
    df_homo = df_homo.sort_values(by ='E-value' )
    
    myCmd_rm = "rm myTmpPeps*"
    os.system(myCmd_rm)
    return (df_homo)

def main():
    args = parse_arguments()
    refArgs = vars(args)

    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFile=refArgs.get("f") ##clusters of peptides
    mydb = refArgs.get("db") ##blast db
    evalue = refArgs.get("evalue") ## evalue cutoff
    
    op = refArgs.get("op") ##Gap opening penalty
    ep = refArgs.get("ep") ##gap extension penalty

    cmd0 = "cut -f 1 " + myFile + " |awk '$0 ~ /^[A-Z]+$/' > myTmp0"
    os.system(cmd0)

    dfTmp = pd.read_csv("myTmp0",header=None,index_col=0,sep="\t",dtype=str)
    df_homo = runBlastp(dfTmp, mydb, evalue, op, ep)
    df_homo = df_homo.rename_axis('Peptides')

    df_homo.to_csv("myBlastp.out.txt", sep="\t", index = True, header=True)
    
    
    #logomatI.to_csv("pssm.txt",sep="\t", index = True, header=True)
    
    cmd1 = "rm myTmp0"
    os.system(cmd1)

if __name__ == '__main__':
    main()



