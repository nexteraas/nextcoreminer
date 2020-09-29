#!/bin/python
import pandas as pd
from venn import venn
import matplotlib.pyplot as plt
import argparse
import sys
import os
import numpy as np

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f',  nargs='*', 
                                    help='The input should be lists of peptides with '
                                    'specificity and frequency ')
   
    my_parser.add_argument('-titles', nargs='*',help='titles of the sets')

    return my_parser.parse_args()


def get_refD(myFiles, myTitles):
   
    refD = {}
    for i in range(0,len(myFiles)):
        files_i = myFiles[i]
        titles_i = myTitles[i]
        cmd0 =  ''' gawk -F "\t" 'length($1)==0 ||$1=="AA_Peptide" || $1 ~ /^[A-Z]+$/' ''' + files_i  + "> myTmp0"  
        os.system(cmd0)
        df_Round = pd.read_csv("myTmp0",header=0,index_col=None,sep="\t")
                
        refD[titles_i] = set(df_Round.iloc[:,0])
    return (refD)


def save_vennPDF(myGroups, refD):
    refMyVenn = {}
    mySize = 0
    for myGroup in myGroups:
        refMyVenn[myGroup] = set(refD.get(myGroup))
        myLen = len(refD.get(myGroup)) 
        if (myLen > 0):
            mySize = myLen
    if (mySize ==0):
        sys.exit("There are no peptides after filterings")
    
    _, top_axs = plt.subplots(ncols=1, nrows=2, figsize=(7, 14))
    
    venn(refMyVenn,  ax =top_axs[0])

    venn(refMyVenn,fmt="{percentage:.1f}%", ax = top_axs[1])
    
    plt.tight_layout()
    plt.savefig("venn.pdf")
    
    plt.close("all")
           

def main():  
    args = parse_arguments()
    refArgs = vars(args)
    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFiles=refArgs.get("f") 
    myTitles = refArgs.get("titles") 
    if (len(myTitles) != len(myFiles)):
        for i in range(0, len(myFiles)):
            myTitles.append("")
    
    if ( len(myFiles)>6 ):
        print ("Venn diagram cannot be shown.")
    elif (len(myFiles)< 2 ):
        sys.exit("Venn diagram needs at least 2 datasets.")
    else:
        refD = get_refD(myFiles, myTitles)
        save_vennPDF(myTitles, refD)


if __name__ == '__main__':
    main()




