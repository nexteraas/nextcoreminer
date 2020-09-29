#!/bin/env python3
import argparse
import pandas as pd
from scipy.spatial import distance
import matplotlib.pyplot as plt
from math import log
import seaborn as sns
import scipy.cluster.hierarchy as hc
import numpy as np
import warnings
warnings.filterwarnings('ignore')


def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='?',  default="test.txt",
                                    help='The input should be a list of peptides with counts, '
                                    'specificity and frequency, by default test.txt is read. ')
    my_parser.add_argument('-l', nargs='?',  default="11",
                                    help='The number of AAs for the peptides, '
                                    'it should be an integer. By default it is 11.')   
    my_parser.add_argument('-o', nargs='?', default="1",
                                    help='Please type an integer if you want to add weight the Bolsum matrix. '
                                    '1: euclidean 2: kullback leibler divergence - relative entropy. '
                                    'By default it is option 1. ') 

    my_parser.add_argument('-colorClass', nargs='?', default="groups",
                                    help='Heatmap color by groups, rounds or cellTypes. ')
    
    return my_parser.parse_args()



def kl_divergence(p, q):
    sum = 0
    for i in range(0,len(p)):    
        sum = sum + (p[i] * log(p[i]/q[i],2) )
    return (sum)


def calcEUDistance (fileTmp, length, blosum_option, AAs, refBlo, colorClass):    
    arrTitle = []    
    for tmp in AAs:
        for leng in range(1,length + 1):
            arrTitle.append("P" + tmp + str(leng))
    refScore = {}   
    myPeps = []
    refScore["Peps"] = ["Classes"] + arrTitle

    file0 = open( fileTmp)  
    while True:
        line0=file0.readline().rstrip()
        if (len(line0) == 0):
            break
        arr = line0.split("\t")
        pep=arr[0]
        if (colorClass == "groups"):
            classes = arr[1]
        elif (colorClass == "rounds"):
            classes = arr[2]
        else:
            classes = arr[3]

        peps = list(pep)
        pep = classes + "_MY_" + pep
        
        if (len(peps) == length):
            score=[]
            for i in range(0,len(peps)):
                aa = peps[i]
                score1 = refBlo.get(aa).split(" ")
                for tmp in score1:
                    tmp2 = ""
                    if (blosum_option == "blosum" ):
                        tmp2 = float(tmp)            ##option1: no frequency ##option3: add freq as a feature ##add weight to each position

                    score.append(str(tmp2)) # ##using the df to calculate EU: 
            refScore[pep] = score 
            myPeps.append(classes + "_MY_" + arr[0] )
    file0.close()

    data={}
    for i in range(0,len(myPeps)):
        for j in range(i,len(myPeps)):
            p=np.array(refScore.get(myPeps[i]),dtype=float)
            q=np.array(refScore.get(myPeps[j]),dtype=float)
            pairs1 = [myPeps[i],myPeps[j]]
            pairs2 = [myPeps[j],myPeps[i]]
            my_dist = distance.euclidean(p, q)
            data[tuple(pairs1)] = my_dist
            data[tuple(pairs2)] = my_dist
    
    dfTmp = pd.DataFrame.from_dict(data, orient='index')

    dfTmp.index = pd.MultiIndex.from_tuples(dfTmp.index.tolist())
    dist_matrix = dfTmp.unstack()
    tuple_list=dist_matrix.columns.tolist()
    first_elements = [a_tuple[1] for a_tuple in tuple_list]
    dist_matrix.columns = first_elements
    dist_matrix =  dist_matrix/max(dist_matrix.max())
    #dist_matrix.to_csv(fileTmp + ".Eudist.txt.my", sep="\t", index = True, header=True)
    return (dist_matrix)
    
    
    
def calcShannonDistance (fileTmp,length, AAs, colorClass):
    file0 = open( fileTmp)  
    myPeps = []
    while True:
        line0=file0.readline().rstrip()
        if (len(line0) == 0):
            break
        arr = line0.split("\t")
        pep=arr[0]
        if (colorClass == "groups"):
            classes = arr[1]
        elif (colorClass == "rounds"):
            classes = arr[2]
        else:
            classes = arr[3]

        peps = list(pep)

        if (len(peps) == length):
            myPeps.append(classes + "_MY_" + arr[0] )
            
    file0.close()

    refBlo={}
    for AA in AAs:
        refBlo[AA] = [0] * length

    for myPep in myPeps:
        peps=list(myPep.split("_MY_")[1])    
        for i in range(0,len(peps)):
            my_arr = refBlo.get(peps[i])
            my_arr[i] = my_arr[i] + 1

    df = pd.DataFrame.from_dict(refBlo, orient="index")/(len(myPeps) + 0.0)

    refScore = {}
    for myPep in myPeps:
        freq = []
        peps=list(myPep.split("_MY_")[1])    
        for i in range(0,len(peps)):
            value = df.at[peps[i],i]
            freq.append(value)

        refScore[myPep] = freq
     
    data={}
    for i in range(0,len(myPeps)):
        for j in range(i,len(myPeps)):
            p=refScore.get(myPeps[i])
            q=refScore.get(myPeps[j])
            rel_entropy = (kl_divergence(p, q) + kl_divergence(q, p))/2
            pairs1 = [myPeps[i],myPeps[j]]
            pairs2 = [myPeps[j],myPeps[i]]
            data[tuple(pairs1)] = rel_entropy
            data[tuple(pairs2)] = rel_entropy
            
    dfTmp = pd.DataFrame.from_dict(data, orient='index')

    dfTmp.index = pd.MultiIndex.from_tuples(dfTmp.index.tolist())
    dist_matrix = dfTmp.unstack()
    tuple_list=dist_matrix.columns.tolist()
    first_elements = [a_tuple[1] for a_tuple in tuple_list]
    dist_matrix.columns = first_elements

    #dist_matrix.to_csv(fileTmp + ".Shannondist.txt", sep="\t", index = True, header=True)
    return (dist_matrix)

def heatmapCluster(df,outName, colorClass, option):
    ##heatmap
    myColors = ["#490092","#24FF24","#FF6DB6", "#009292", "#396AB1","#DA7C30", "#3E9651", "#CC2529", "#535154", "#6B4C9A", "#922428", "#948B3D", "#984ea3", "#ff7f00", "#ffff33", "#a65628","#e1f5fe","#1A4756","#33918D","#F6D82D","#E39517","#D4342D"]
    myColumns = df.columns
    myRows = df.index
    myColumns = [i.split('_MY_')[0] for i in myColumns] 
    myRows = [i.split('_MY_')[0] for i in myRows] 
    df.columns = myColumns
    df.index = myRows
    mycols = []
    refCol={}
    for column in myColumns:
        if (column not in refCol):
            refCol[column] = myColors[0]
            myColors.remove(myColors[0])
            
    colors=pd.Series(df.columns,index=df.columns).map(refCol)
    colors1=pd.Series(df.columns).map(refCol)

    colors1=colors1.rename(colorClass)
    
    
    df.columns = pd.RangeIndex(df.columns.size)
        
    linkage = hc.linkage(df, method="complete")
    
    #check cmap='BuPu,YlGnBu, GnPu, Blues    vlag'  
    fig = plt.figure()
    myTitle="Kullback–Leibler Divergence"
    if (option=="1"):
        myTitle="Euclidean Distance"
    

    g = sns.clustermap(df, row_linkage=linkage, col_linkage=linkage, xticklabels=False, yticklabels=False, linecolor='gray',col_colors=colors1,cmap="BuPu",figsize=(8,8))
    g.fig.suptitle(myTitle)

    for (label, colour) in refCol.items():
        g.ax_row_dendrogram.bar(0.5, 0.0, color=colour, label="{}".format(label))

    nrows = 25
    ncols = int(np.ceil(len(refCol) / float(nrows)))
    legend = g.ax_row_dendrogram.legend(ncol=ncols, title=colorClass,bbox_to_anchor=(5.3, 1), loc='upper left',fontsize='x-small')
    legend.get_frame().set_facecolor('none')
    #legend.get_frame().set_linewidth(0.0)
    
    #g.cax.set_position([0.92, .29, .03, .4])
    plt.savefig(outName + ".png",dpi=600,bbox_inches='tight')
    plt.close('all')
    
def read_blossum_matrix(fileBloTmp):
    AA = ""
    fileBlo = open(fileBloTmp)
    refBlo = {}
    index=20
    
    blosumAA = list("ARNDCQEGHILKMFPSTWYV")
    while True:
        line = fileBlo.readline().rstrip()
        if (len(line) == 0):
            break
        if ("#" not in line):
            arr = line.split()
            countAA = 0
            for i in range(0,len(arr)):
                if (arr[i] in blosumAA):
                    countAA = countAA + 1
            if (countAA > 1):
                for i in range(0,len(arr)):
                    if (arr[i]  in "BZX*"):
                        index = i
                        break
                    else:
                        AA = AA + arr[i]

            else:
                if (arr[0] not in "BZX*"):
                    score1 = []
                    for i in range(1,index + 1):                    
                        score = arr[i]                    
                        score1.append(score)
                    refBlo[arr[0]] = " ".join(score1)
    fileBlo.close()

    AAs = list(AA)
    return (refBlo, AAs)

def main():
    blosum_option="blosum"
    
    refBlo, AAs = read_blossum_matrix("/opt/nextera/galaxy/data/custom/BLOSUM62.txt")

    args = parse_arguments()
    refArgs = vars(args)

    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFile = refArgs.get("f") ##peptide file with frequency
    myLength = int(refArgs.get("l")) ##peptide sequence length
    colorClass = refArgs.get("colorClass") ##Heatmap color by classes, rounds or cellTypes
    option = refArgs.get("o") #option for euclidean distance or relative entropy

    if (option == "1"):
        heatmapCluster(calcEUDistance(myFile, myLength, blosum_option, AAs, refBlo, colorClass),"HeatmapWithinCluster", colorClass, option)
    else:
        heatmapCluster(calcShannonDistance(myFile, myLength, AAs, colorClass),"HeatmapWithinCluster", colorClass, option)

if __name__ == "__main__":
    main()

