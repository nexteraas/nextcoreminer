#!/bin/python

import argparse
import seaborn as sns
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt
from math import ceil
import os
import sys
import numpy as np
import re

sns.set(font_scale=2)
plt.rcParams["font.family"] = "DejaVu Sans"

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f',  nargs='*', 
                                    help='The input should be lists of peptides ')
   
    my_parser.add_argument('-colors', nargs='*',help='colorbars of the plots')
    my_parser.add_argument('-titles', nargs='*',help='titles of the plots')
    my_parser.add_argument('-pos',  nargs='?', 
                                    help='The input should be lists of AA positions. '
                                    'If no -p option, then all the AAs are included.')
    my_parser.add_argument('-optionScale', nargs='?', default="no",help='Do you want to scale AA prevalences by read counts?')
    my_parser.add_argument('-optionOut', nargs='?', default="no",help='Do you want to ignore some positions?')

    return my_parser.parse_args()

def rd(x,y=0):
    ''' A classical mathematical rounding by Voznica '''
    m = int('1'+'0'*y) # multiplier - how many positions to the right
    q = x*m # shift to the right by multiplier
    c = int(q) # new number
    i = int( (q-c)*10 ) # indicator number on the right
    if i > 0:
        c += 1
    return c/m


def getLogomatI (df_Round, pos):
    logomatI = lm.alignment_to_matrix(sequences=list(df_Round.iloc[:,0]))
    
    ## 1st position should not be 1
    if (min(list(logomatI.index)) == 0):
        logomatI.index = logomatI.index + 1

    if (pos is not None):
        positions = [int(i) for i in pos]
        for pos in positions:
            logomatI.loc[pos]=0
    
    logomatCount = logomatI
    logomatI = logomatI.transpose()/df_Round.shape[0]
    
    logomatI.index.name = None
    
    maxForRow = logomatI.max(axis=0)
    globalMax = maxForRow.max()
    
    maxVal = rd(globalMax,1)
    
    return (logomatCount, logomatI, maxVal)
    


def plot_Heatmap(df, pos, g, index, ncols,shrink_value, myTitle, myColor, maxValue):
    logomatCountTmp, logomat_Round, maxVal = getLogomatI(df, pos)
    
    autosize = True
    if (logomat_Round.shape[1] == 11):
        shrink_value = shrink_value - 0.13
    g = sns.heatmap(logomat_Round, cmap = myColor,  vmin= 0, vmax=maxValue,square=True,linewidth=0.5,cbar_kws = dict(orientation='horizontal',pad=0.01,shrink=shrink_value, label= 'Proportion of sequences'))
    g.tick_params(axis='both', which='both', length=0)
    g.set_title(myTitle)

    #bottom, top = g.get_ylim()
    #g.set_ylim(bottom + 0.5, top - 0.5)
    for _, spine in g.spines.items():
        spine.set_visible(True)

    g.xaxis.tick_top()
    g.xaxis.set_label_position('top')
    
    first_column = ( (index -1) % ncols == 0 )
    if first_column:
        g.set_ylabel("Amino acid")
        
    first_row = ( index <= ncols  )    
    if first_row:
        g.set_xlabel("Peptide position")
        
    for tick in g.get_yticklabels():
        tick.set_rotation(0)    
    

def paramReturn (_max): 
    
    ncols = 3
    if (_max == 1 ):
        f = plt.figure(figsize=(10,12))
        shrink_value = 0.55
        ncols = _max
    elif (_max == 2):
        f = plt.figure(figsize=(14,12))
        shrink_value = 0.85
        ncols = _max
    elif (_max == 3):
        f = plt.figure(figsize=(20,12))
        ncols = _max
        shrink_value = 0.92
    elif (_max == 4) :
        f = plt.figure(figsize=(30,32))
        shrink_value = 0.81
    elif (_max == 5):
        f = plt.figure(figsize=(30,32))
        shrink_value = 0.78
    elif (_max == 6):
        f = plt.figure(figsize=(30,32))
        shrink_value = 0.75
    elif (_max == 7):
        f = plt.figure(figsize=(30,52))
        shrink_value = 0.94
    elif (_max == 8):
        f = plt.figure(figsize=(30,52))
        shrink_value = 0.9
    else:
        f = plt.figure(figsize=(30,52))
        shrink_value = 0.88

    return f, ncols, shrink_value


def get_dfRound(files_i, optionScale):
    cmd0 =  ''' gawk -F "\t" 'length($1)==0 ||$1=="AA_Peptide" || $1 ~ /^[A-Z]+$/' ''' + files_i  + "> myTmp0"  
    os.system(cmd0)
    df_Round = pd.read_csv("myTmp0",header=0,index_col=None,sep="\t")
    if (optionScale == "yes"):
        all_columns = list(df_Round.columns)
        col_count = None
        for col in all_columns:
            if ("Count" in col ):
                col_count = col
        if (col_count is None):
            sys.exit("You have chosen to scale aa prevalences by read counts. However, there is no count column in your file.")
        else:

            newdf = pd.DataFrame(np.repeat(df_Round.values,df_Round[col_count].astype(int),axis=0))
            newdf.columns = df_Round.columns
            df_Round = newdf
    return (df_Round)

def main():
    
    args = parse_arguments()
    refArgs = vars(args)
    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFiles=refArgs.get("f") ##rounds file
    if (len(myFiles) >9 ):
        sys.exit("There should be max. 9 files.")

    myTitles = refArgs.get("titles") ## titles of the plots

    if (len(myTitles) != len(myFiles)):
        for i in range(0, len(myFiles)):
            myTitles.append("")
            
    myColors = refArgs.get("colors") ## colors of the plots
    if (len(myColors) != len(myFiles)):
        for i in range(0, len(myFiles)):
            myColors.append("Reds")
    
    
    optionScale = refArgs.get("optionScale") ## scale AA prevalence with counts
    optionOut = refArgs.get("optionOut") ## Do you want to ignore some positions
    
    pos = None
    if (optionOut == "yes"):
        if (refArgs.get("pos") is None):
            sys.exit("You have chosen to ignore some positions. However, there are no positions provided")
        else:
            position = refArgs.get("pos") ##positions to be sank out
    
            delimiters = ["X",",", ";", " "]
            regexPattern = '|'.join(map(re.escape, delimiters))
            positionTmp= re.split(regexPattern, position)
            pos = [x for x in positionTmp if x] #position for each position
    
    _max = len(myFiles)
    fig, ncols, shrink_value = paramReturn(_max)
    nrows = ceil(_max/ncols)
    
    ##To get the max value of the matrix
    ##the code below was written twice
    arrValue = []
    for i in range(0,len(myFiles)):
        df_Round = get_dfRound(myFiles[i], optionScale)
        
        logomatCount, dfTmp, maxValue = getLogomatI (df_Round, pos)
        arrValue.append(maxValue)
        
        myTitle = myTitles[i]
        logomatCount.index.name = myTitle
        logomatCount.to_csv("AAPrevalenceCounts.txt",sep="\t", index = True, header=True,mode='a')
        
    maxValue = max(arrValue)
    
    
    for i in range(0,len(myFiles)):
        df_Round = get_dfRound(myFiles[i], optionScale)
        
        myTitle = myTitles[i]
        myColor = myColors[i]
        axis = fig.add_subplot(nrows,ncols,i+1)
        plot_Heatmap(df_Round, pos, axis, (i + 1),  ncols, shrink_value, myTitle, myColor, maxValue)
        
        cmd1 = "rm myTmp0"
        os.system(cmd1)
        
    fig.tight_layout()
    plt.yticks(rotation=0) 
    plt.savefig( "AAPrevalencePlot.pdf")
    plt.close('all')

           
if __name__ == '__main__':
    main()





