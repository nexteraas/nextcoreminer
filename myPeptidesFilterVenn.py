#!/bin/python
import itertools
import pandas as pd
from venn import venn
import matplotlib.pyplot as plt
import argparse
import sys
import os
from math import ceil
import numpy as np
import re
from myHomoPeps import getHomoPeptides


def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f',  nargs='*', 
                                    help='The input should be lists of peptides with counts and frequency ')
   
    my_parser.add_argument('-g', nargs='?',help='round_exp groups for venn')
    my_parser.add_argument('-filterFreq', nargs='?',default="no", help='remove peptides if R1_freq < R2_freq')
    my_parser.add_argument('-filterHomo', nargs='?',default="no", help='remove peptides if they do not share homology with epitope')
    my_parser.add_argument('-filterNegCells', nargs='?',default="no", help='remove peptides if they appear in the negative panning cells')
    my_parser.add_argument('-negCells', nargs='?',default="",help='negative cells, e.g. R1_380 R2_380')
    my_parser.add_argument('-negCellsCounts', nargs='?',default="0", help='remove peptides if they appear in the negative panning cells with counts higher than * number')
    my_parser.add_argument('-epi',  nargs='?', default="KLLTQHFVQENY",
                                    help='The input should be the epitope of interest.')
    my_parser.add_argument('-o', nargs='?', default="Identity",help='Identity or Similarity to measure homology  ')

    my_parser.add_argument('-threshold', nargs='?',  default="80",
                                    help='The cutoff of identity/similarity for fasta.')
    my_parser.add_argument('-n', nargs='?',  default="1",
                                    help='The cutoff of number of overlap AA')
    my_parser.add_argument('-op', nargs='?',  default="100", help='Gap opening penalty')
    my_parser.add_argument('-ep', nargs='?',  default="0.5", help='Gap extension penalty')
    my_parser.add_argument('-c', nargs='?',  default="0",
                                    help='Only include the peptides with counts higher than a number. '
                                    'By default it is 0, meaning all peptides will be included')
    my_parser.add_argument('-filterCounts', nargs='?',default="no", help='remove peptides from neg cells with counts higher than XX')
    
    return my_parser.parse_args()


def get_petal_labels(datasets, fmt="{size}"):
    """Generate petal descriptions for venn diagram based on set sizes"""
    datasets = list(datasets)
    n_sets = len(datasets)
    dataset_union = set.union(*datasets)
    universe_size = len(dataset_union)
    petal_labels = {}
    
    refSet = {}
    
    for logic in get_logics(n_sets):
        included_sets = [
            datasets[i] for i in range(n_sets) if logic[i] == "1"
        ]
        excluded_sets = [
            datasets[i] for i in range(n_sets) if logic[i] == "0"
        ]
        petal_set = (
            (dataset_union & set.intersection(*included_sets)) -
            set.union(set(), *excluded_sets)
        )
        petal_labels[logic] = fmt.format(
            logic=logic, size=len(petal_set),
            percentage=(100*len(petal_set)/max(universe_size, 1))
        )
        
        
        if (len(petal_set) >0):
            refSet[logic] = list(petal_set)
    
    return refSet

def get_logics(n_sets):
    """Generate intersection identifiers in binary (0010 etc)"""
    for i in range(1, 2**n_sets):
        yield bin(i)[2:].zfill(n_sets)


def checkValueIncrease(arr_tmp): 
    
    res = all(i <= j for i, j in zip(arr_tmp, arr_tmp[1:])) 

    return (res)
    
    
def filtering_Epi(files,filterCounts, counts, myGroups, filterHomo, myEpi, option, threshold, num_overlap, filterNegCells, negCells, negCellsCounts, op, ep):
    my_df = pd.DataFrame()

    refLog = {}
    d = {}
    for myFile in files:
        df0 =pd.read_csv(myFile,header=None)
        name2 = df0.iloc[1].apply(str).values[0].replace("#","")
        d[name2] = pd.read_csv(myFile,index_col=0, header=0,sep="\t",dtype=str,skiprows=[0,1,2,3,4,5])
        d[name2].columns = [name2 + "_" + str(col).replace("_","").replace(" ","")  for col in d[name2].columns ]
        
        
        if (name2 in myGroups):
            refLog["Original" + "AND" + name2] = d[name2].shape[0]
        
        ######### START of Filtering4 - Counts higher than a number   #######
        if (filterCounts == "yes"):
            cond_count = d[name2][name2 + "_TotalCount"].astype(float) > counts
            d[name2] = d[name2][cond_count]

            if (name2 in myGroups):
                refLog["filterCounts" + "AND" + name2] = d[name2].shape[0]
        
        ######## END of Filtering 4  #######
        
        
        ######### START of Filtering 1 - Homology                             #######
        ######### Filter 1: retrive peptides sharing homology with eptiope    #######
        if (filterHomo == "yes"):
            my_df_homo = getHomoPeptides(d[name2], myEpi,option, threshold, num_overlap, op, ep)
            d[name2] = d[name2][d[name2].index.isin(my_df_homo.index)]
            
            cmd_remove = "rm myFasta_out.txt" 
            os.system(cmd_remove)
            
            if (name2 in myGroups):
                refLog["filterHomo" + "AND" + name2] = d[name2].shape[0]
            
        my_df = pd.concat([my_df,d[name2]],axis=1,sort=True)
        ######### END of Filtering 1  #######

    ######### START of Filtering3 - Negative cells   #######
    if (filterNegCells == "yes"):
        
        negCellsCounts_new = ['%s_TotalCount' % (x) for x in negCells]
        
        condition1 = ~ (my_df[negCellsCounts_new].isna().all(axis = 1))  ## should not be NA in ALL the neg cells
        condition2 = (my_df[negCellsCounts_new].astype(float) > negCellsCounts).any(axis = 1)  ## if in negative cells, then remove those with counts higher than a number

        dfTmp = my_df[condition1 & condition2]
        
        my_df = my_df[~ (my_df.index.isin (list(dfTmp.index)) )]    
        my_df = my_df[3:].dropna(how='all',axis=0)
    
        
        for myGroup in myGroups:
            refLog["filterNeg" + "AND" + myGroup] = my_df[my_df[myGroup + "_Frequency"].notnull()].shape[0]
        
    ######### END of Filtering3    ########
    
    return my_df, d, refLog


def get_rounds_exps(files,d):
    refD={}
    for myFile in files:
        df0 = pd.read_csv(myFile,header=None)
        name2 = df0.iloc[1].apply(str).values[0].replace("#","")
        exec( "%s = frozenset(d.get(name2).dropna().index.tolist())" % (name2  ))
        txt = "%s" % (name2)

        refD[txt] = set(vars()[txt])

    refExps={}
    refRounds={}
    for keys in refD:
        if ("R0_Amp" not in keys):
            exp=keys.split("_")[1]
            rou=keys.split("_")[0]
            refExps[exp] = 1
            refRounds[rou]=1

    rounds = sorted(refRounds)
    exps   = sorted(refExps)
    
    return (refD, rounds, refExps, exps)


def merge_df(refExps, rounds, my_df,filterFreq, refD, myGroups):
    df_merged = pd.DataFrame()

    for exp in refExps:
        
        arr = ["R0_Amp_Rank","R0_Amp_Frequency","R0_Amp_TotalCount"]
        arr_exp = []
    
        for rou in rounds:
            comb = rou + "_" + exp
            comb1 = comb + "_" + "Rank"
            comb2 = comb + "_" + "Frequency"
            comb3 = comb + "_" + "TotalCount"
            arr = arr + [comb1,comb2,comb3]
            arr_exp = arr_exp + [comb1,comb2,comb3]
        
        df = my_df[arr]
        freq=[]
        count=[]
        for col in list(df.columns):
            if ("freq" in col.lower()):
                freq.append(col)
            if ("count" in col.lower()):
                count.append(col)
                
        
        ######### START of Filtering 2 #######
        ######### Filter 2: remove peptides with in lower rounds with higher frequency. i.e. R1_freq should be < R2_freq  
       
        if (filterFreq == "yes"):
            
            cond1 = df[freq].fillna("-1").apply(lambda row : checkValueIncrease(np.array(list(row),dtype=float)), axis = 1)  ##peptides should not be NA in all non-R0, and freq should increase 
            
            df = df[cond1  ]
           
        ######### END of Filter 2 #########  
        
        
        df = df[arr_exp]
        
        ##filter those with NA for all the rows
        df = df.dropna(how='all',axis=0)
                
        df_merged = pd.concat([df_merged,df],axis=1,sort=True)
        
        for rou in rounds:
            myPeps = set(list(df[df[rou + "_" + exp + "_Frequency"].notnull()].index))
            refD[rou + "_" + exp] = myPeps
    
    refLog = {}  
    if (filterFreq == "yes"):
        for myGroup in myGroups:
            refLog["filterFreq" + "AND" + myGroup] = df_merged[df_merged[myGroup + "_Frequency"].notnull()].shape[0]
            
    return df_merged, refD, refLog

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
    if (len(myGroups)<2 ):
        print("Venn diagram needs at least 2 datasets.")
    elif(len(myGroups)>6):
        print("Venn diagram cannot be shown for more than 6 datasets.")
    else:
        _, top_axs = plt.subplots(ncols=1, nrows=2, figsize=(7, 14))
        venn(refMyVenn,  ax =top_axs[0])
        venn(refMyVenn,fmt="{percentage:.1f}%", ax = top_axs[1])

        plt.tight_layout()
        plt.savefig("venn.pdf")
        plt.close("all")

    
    ##OUTPUT THE VENN GROUPS
    out = open("venn.txt","w")
    out.write("peptide" + "\t" + "group" + "\t" + "rounds" + "\t" + "cells" + "\n")
    refSet = get_petal_labels(refMyVenn.values())
    tmpGroups = list(refMyVenn)
    
    for keys in refSet:
        logics = list(keys)
        nameGroups = []
        for i in range(0,len(logics)):
            if (logics[i] == "1"):
                nameGroups.append(tmpGroups[i] )
        
        for tmp in refSet.get(keys):
            
            refRounds={}
            refCells ={}
            for tmp_t in nameGroups:
                arr_t = tmp_t.split("_")
                refRounds[arr_t[0]] = 1
                refCells[arr_t[1]] = 1
            rounds = sorted(refRounds)
            cells  = sorted(refCells)
            
            out.write(tmp + "\t" + "+".join(nameGroups) + "\t" + "+".join(rounds) + "\t" + "+".join(cells) + "\n")
    out.close()
       

def main():  
    args = parse_arguments()
    refArgs = vars(args)
    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    files=refArgs.get("f") 
    myGroup = refArgs.get("g") 
    filterFreq = refArgs.get("filterFreq") 
    filterHomo = refArgs.get("filterHomo") 
    filterNegCells = refArgs.get("filterNegCells") 
    filterCounts = refArgs.get("filterCounts") 
    negCell = refArgs.get("negCells") 
   
    negCells = []
    delimiters = ["X",",", ";", " "]
    regexPattern = '|'.join(map(re.escape, delimiters))
    if (len(negCell) != 0):
        negCellTmp= re.split(regexPattern, negCell)
        negCells = [x for x in negCellTmp if x] 
    
    myGroupTmp= re.split(regexPattern, myGroup)
    myGroups = [x for x in myGroupTmp if x] 
    
    negCellsCounts = int(refArgs.get("negCellsCounts"))
    
    myEpi = refArgs.get("epi") 
    option = refArgs.get("o") ##option of identity or threshold
    threshold = float(refArgs.get("threshold")) ## threshold cutoff
    num_overlap = float(refArgs.get("n")) ## number of overlap AA
    
    op = refArgs.get("op") ##Gap opening penalty
    ep = refArgs.get("ep") ##gap extension penalty
    counts = int(refArgs.get("c")) ## if it is 0, then all peptides are included
    
    
    ## functions start:
    my_df,d, refLog1 = filtering_Epi(files, filterCounts, counts, myGroups, filterHomo, myEpi,option, threshold, num_overlap, filterNegCells, negCells, negCellsCounts, op, ep)

    refD, rounds, refExps, exps = get_rounds_exps(files, d) 

    df_merged, refD, refLog2 = merge_df(refExps, rounds,my_df,filterFreq, refD, myGroups )

    ## output log
    outLog = open("FilterLog.txt","w")
    outLog.write("python3 myPeptidesFilterVenn.py ")
    n = len(sys.argv) 
    
    for i in range(1, n): 
        outLog.write(sys.argv[i] + " ")
    outLog.write("\n" + "\n" )
    
    outLog.write("Filter" + "\t" + "Group" + "\t" + "#PeptidesRemained" + "\n" )
    for myGroup in myGroups:
        for keys in refLog1:
            if (keys.split("AND")[1] == myGroup):
                outLog.write("\t".join(keys.split("AND")) + "\t" + str(refLog1.get(keys)) + "\n")
        for keys in refLog2:
            if (keys.split("AND")[1] == myGroup):
                outLog.write("\t".join(keys.split("AND")) + "\t" + str(refLog2.get(keys)) + "\n")
    
    outLog.close()


    save_vennPDF(myGroups, refD)


    arr_Groups = []
    for myGroup in myGroups:
        comb1 = myGroup + "_" + "Rank"
        comb2 = myGroup + "_" + "Frequency"
        comb3 = myGroup + "_" + "TotalCount"
        arr_Groups.append(comb1)
        arr_Groups.append(comb2)
        arr_Groups.append(comb3)

    df = df_merged[arr_Groups]
    df = df.dropna(how='all',axis=0)  
    df2 = pd.read_csv("venn.txt",index_col = 0, header = 0,sep="\t",dtype = str)
    new_df  = pd.concat([df2,df],axis = 1,sort = True)
    df = new_df

    freq=[]
    count=[]
    for col in list(df.columns):
        if ("freq" in col.lower()):
            freq.append(col)
        if("count" in col.lower()):
            count.append(col)

    ##take the max value if peptides appear in > dataset
    df['MaxFreq']  = df[freq].astype(float).max(axis=1)
    df['MaxCount'] = df[count].astype(float).max(axis=1)        

    df.fillna("NA").to_csv("peptide_matrix.txt", sep="\t", index = True, header=True)
        

if __name__ == '__main__':
    main()






