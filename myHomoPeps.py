#!/bin/env python3
import argparse
import os
import sys
import pandas as pd
import numpy as np
from math import ceil


def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='?',  help='a file containing peptides on the 1st column ')

    my_parser.add_argument('-epi',  nargs='?', default="",
                                    help='The input should be the epitope of interest.')

    my_parser.add_argument('-o', nargs='?', default="Identity",help='Please select Identity or Similarity  ')
    my_parser.add_argument('-threshold', nargs='?',  default="80",
                                    help='The cutoff of identity/similarity for fasta. By default it is 80.')
    my_parser.add_argument('-op', nargs='?',  default="100", help='Gap opening penalty')
    my_parser.add_argument('-ep', nargs='?',  default="0.5", help='Gap extension penalty')

    my_parser.add_argument('-n', nargs='?',  default="1",
                                    help='The cutoff of number of overlap AA')

    return my_parser.parse_args()


def getHomoPeptides (dfTmp, myEpi, option, threshold, num_overlap, op, ep):
    df_homo =  pd.DataFrame()
    myCmd1 = '''echo ">" '''+ myEpi + '''"\n"''' + myEpi + "> myDatabaseEpi.fa"
    os.system(myCmd1)
        
    fasta36 = "/usr/bin/fasta-36.3.8h/bin/fasta36"

    maxi = 60000
    num_files = ceil(len(list(dfTmp.index))/maxi)

    arrTmps = np.array_split(list(dfTmp.index), num_files)   
    for i in range(0,len(arrTmps)):
        tmps = arrTmps[i]
        outTmp = open("myDatabase" + str(i) + ".fa", "w")
        for tmp in tmps:
            outTmp.write(">" + tmp + "\n")
            outTmp.write(tmp + "\n")
        outTmp.close()

        myCmd2 = fasta36 + " -s BP62 myDatabaseEpi.fa myDatabase"  + str(i) + ".fa -f " + op + " -g " + ep  + " -b=" + str(maxi) + " > myDatabase" + str(i) + ".out"
        os.system(myCmd2)
        
        myCmd3 = '''gawk 'length($0)>0' myDatabase''' + str(i) + ".out > myDatabase.tmp"
        os.system(myCmd3)
        
        file0 = open("myDatabase.tmp" )
        arrLine = []
        while True:
            line = file0.readline().rstrip()
            if (len(line) == 0):
                break
            arrLine.append(line)
        file0.close()

        out = open("myDatabase" + str(i) + ".final","w")
        for index in range(0,len(arrLine)):
            if(arrLine[index].startswith(">>")):
                pep = arrLine[index].replace(">>","").split()[0]
                identity = ""
                similarity = ""
                num = 0
                num_db = 0
                num_sd = 0
                identity = arrLine[index + 2].split("; ")[1].split(" identity (")[0].lstrip().replace("%","") 
                similarity = arrLine[index + 2].split("; ")[1].split(" identity (")[1].split(" similar")[0].lstrip().replace("%","")
                dotsLine = arrLine[index + 5].lstrip()
                num = len(dotsLine)
                num_db = dotsLine.count(":")
                num_sd = dotsLine.count(".") + num_db

                out.write(pep + "\t" + identity + "\t" + similarity + "\t" + str(num) + "\t" + str(num_db) + "\t" + str(num_sd) + "\n")

        out.close()      

        df_tmp = pd.read_csv("myDatabase" + str(i) + ".final", index_col=0, header=None,sep="\t")
        df_homo = pd.concat([df_homo, df_tmp],axis=0,sort=True)
    
    myCmd_cat = "cat myDatabase*.out > myFasta_out.txt"
    os.system(myCmd_cat)
    myCmd6 = "rm myDatabase*" 
    os.system(myCmd6)
    columnNames = [ "Identity(%)" , "Similarity(%)" , "#OverlapAA" ,"#IdentityAA" , "#SimilarityAA" ]
    df_homo.columns = columnNames
    
    cond1 = df_homo[option + "(%)"] >= threshold
    cond2 = df_homo["#OverlapAA"] >= num_overlap

    df_homo = df_homo[cond1 & cond2 ]
    return (df_homo)

def main():
    args = parse_arguments()
    refArgs = vars(args)

    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFile=refArgs.get("f") 
    myEpi = refArgs.get("epi") ##epitope of interest
    option = refArgs.get("o") ##option of identity or threshold
    threshold = float(refArgs.get("threshold")) ## threshold cutoff
    num_overlap = float(refArgs.get("n")) ## number of overlap AA
    op = refArgs.get("op") ##Gap opening penalty
    ep = refArgs.get("ep") ##gap extension penalty

    cmd0 = "cut -f 1 " + myFile + " |awk '$0 ~ /^[A-Z]+$/' > myTmp0"
    os.system(cmd0)

    dfTmp = pd.read_csv("myTmp0",header=None,index_col=0,sep="\t",dtype=str)
    df_homo = getHomoPeptides(dfTmp, myEpi, option, threshold, num_overlap,op,ep)
    df_homo = df_homo.rename_axis('Peptides')

    df_homo.to_csv("myHomoPeps.txt", sep="\t", index = True, header=True)
    cmd1 = "rm myTmp0"
    os.system(cmd1)

if __name__ == '__main__':
    main()


