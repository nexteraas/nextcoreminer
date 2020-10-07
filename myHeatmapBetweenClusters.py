#!/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as hc
import logomaker as lm
import math
import numpy as np
import warnings
warnings.filterwarnings('ignore')


def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='*',  
                                    help='The input should be a list of peptides on the 1st column  ')
        
    return my_parser.parse_args()


def calcEUDistance (myFiles):    
    AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
       'R', 'S', 'T', 'V', 'W', 'Y']

    refClusters = {}
    for myFile in myFiles:
        df = pd.read_csv(myFile,sep="\t",header=0,index_col=0)
        pep = list(df.index)

        logomat  = lm.alignment_to_matrix(sequences=pep)
        
        ###### fill in 0 for the missing AA, so that there are alaways 20 AAs
        columnAAs = logomat.columns
        for AA in AAs:
             if (AA not in columnAAs):
                logomat[AA] = 0

        columnAAs_sorted = sorted(list(logomat.columns) )
        logomat = logomat[columnAAs_sorted]
        
        logomatI = lm.transform_matrix(df = logomat, from_type='counts', to_type='information')
        
        fileName = (df.index.name).split("_")[0]
        refClusters[fileName] = logomatI

    data={}
    keys = list(refClusters.keys())
    
    for i in range (0,len(keys)):
        for j in range(i, len(keys)):
            df1 = refClusters.get(keys[i]) 
            df2 = refClusters.get(keys[j]) 
            
            df3 = df1 - df2
            df3_squre = df3**2
            sumAll = df3_squre.values.sum()
            my_dist = math.sqrt(sumAll)
           
            pairs1 = [keys[i],keys[j]]
            pairs2 = [keys[j],keys[i]]

            data[tuple(pairs1)] = my_dist
            data[tuple(pairs2)] = my_dist

    dfTmp = pd.DataFrame.from_dict(data, orient='index')
    dfTmp.index = pd.MultiIndex.from_tuples(dfTmp.index.tolist())
    dist_matrix = dfTmp.unstack()
    tuple_list=dist_matrix.columns.tolist()
    first_elements = [a_tuple[1] for a_tuple in tuple_list]
    dist_matrix.columns = first_elements
    dist_matrix =  dist_matrix/max(dist_matrix.max())
    
    dist_matrix.to_csv("Distance_matrix.txt", sep="\t", index = True, header=True)
    return (dist_matrix)
    
    
#####the same as myHeatMapWithinCluster.py
def heatmapCluster(df,outName, colorClass):
    ##heatmap
    myColors = ["#490092","#24FF24","#FF6DB6", "#009292", "#396AB1","#DA7C30", "#3E9651", "#CC2529", "#535154", "#6B4C9A", "#922428", "#948B3D", "#984ea3", "#ff7f00", "#ffff33", "#a65628","#e1f5fe","#1A4756","#33918D","#F6D82D","#E39517","#D4342D"]
    myColumns = df.columns
    myRows = df.index
    
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
    
    g = sns.clustermap(df, row_linkage=linkage, col_linkage=linkage, xticklabels=False, yticklabels=False, linecolor='gray',col_colors=colors1,cmap="BuPu",figsize=(8,8))
    #g.fig.suptitle(myTitle)

    for (label, colour) in refCol.items():
        g.ax_row_dendrogram.bar(0.5, 0.0, color=colour, label="{}".format(label))

    nrows = 25
    ncols = int(np.ceil(len(refCol) / float(nrows)))
    legend = g.ax_row_dendrogram.legend(ncol=ncols, title=colorClass,bbox_to_anchor=(5.3, 1), borderaxespad=1., loc='upper left',fontsize='x-small')
    legend.get_frame().set_facecolor('none')
    
    plt.savefig(outName + ".png",dpi=600,bbox_inches='tight')
    plt.close('all')
    

def main():
    
    args = parse_arguments()
    refArgs = vars(args)

    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFiles = refArgs.get("f") ##peptide file with frequency
    
    colorClass = "Clusters"
    heatmapCluster(calcEUDistance(myFiles),"HeatmapBetweenClusters", colorClass)
    
if __name__ == "__main__":
    main()



