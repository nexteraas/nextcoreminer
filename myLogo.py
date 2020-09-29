#!/bin/python
import pandas as pd
import matplotlib.pyplot as plt
import logomaker as lm
import argparse
import os
import sys
import re

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='*', 
                                    help='The input should be a list of peptides with counts, '
                                    'specificity and frequency. By default test.txt is read. ')
    my_parser.add_argument('-pos',  nargs='?', 
                                    help='The input should be lists of AA positions. '
                                    'If no -p option, then all the AAs are included.')
    my_parser.add_argument('-titles',  nargs='*',default="", 
                                    help='The title of the logo plot. ')
    
    my_parser.add_argument('-type', nargs='?',  default="logo",
                                    help="logo or ice-logo")
    my_parser.add_argument('-optionInput', nargs='?',  default="peptide_list", help='input file option')

    my_parser.add_argument('-optionOut', nargs='?', default="no",help='Do you want to ignore some positions?')

    return my_parser.parse_args()


def showlogo(df,myTitle, myLogoType, axis):
    mylogo=lm.Logo(df,font_name="DejaVu Sans", ax=axis)
    mylogo.style_spines(visible=False)
    mylogo.style_spines(spines=["left","bottom"],visible=True)
    if (myLogoType == "logo"):
        mylogo.ax.set_ylabel("Information (bits)")
    else:
        mylogo.ax.set_ylabel("Weight (bits)")

    mylogo.ax.set_title(myTitle)
    mylogo.ax.xaxis.set_ticks_position("none")
    mylogo.ax.xaxis.set_tick_params(pad=-1)
    mylogo.ax.grid(False)
    mylogo.fig.tight_layout()

def readPSSM(pssm):
    AAs = list("ARNDCQEGHILKMFPSTWYV")

    out = open("myTmp.PSSM.matrix","w")
    with open(pssm) as f_in:   
        for line0 in f_in:
            line = line0.rstrip()
            if len(line.split()) == 0:
                continue
            else:
                if ( not (line.startswith("Last position-specific scoring"))):
                    arrLine = line.split()
                    #print (arr)
                    if ("Lambda" not in arrLine and "Gapped" not in arrLine and "Ungapped" not in arrLine):
                        if (arrLine[0] in AAs):
                            out.write("\t" + "\t".join(arrLine[0:20]) + "\n")
                        else:
                            out.write(arrLine[0] + "\t" + "\t".join(arrLine[2:20 + 2]) + "\n")
    out.close()                

def getLogomatI (optionInput, pos, myLogoType, myFile):
    if (optionInput == "peptide_list"):
    
        cmd0 = "cut -f 1 " + myFile + " |awk '$0 ~ /^[A-Z]+$/' > myTmp0"
        os.system(cmd0)
        peps = pd.read_csv("myTmp0",index_col=False, header=None,sep="\t",dtype=str)

        logomat  = lm.alignment_to_matrix(sequences=peps.iloc[:,0])        

        if (myLogoType == "logo"):
            logomatI = lm.transform_matrix(df = logomat, from_type='counts', to_type='information')
        else:
            logomatI = lm.transform_matrix(df = logomat, from_type='counts', to_type='weight')
            logomatI = lm.transform_matrix(df = logomatI,center_values=True)
        cmd1 = "rm myTmp0 "
        os.system(cmd1)
        
    else:
        readPSSM(myFile)
        myLogoType = "ice-logo"
        logomatI = pd.read_csv("myTmp.PSSM.matrix",index_col=0, header=0,sep="\t")
        logomatI = lm.transform_matrix(df = logomatI,center_values=True)
        cmd2 = "rm myTmp.PSSM.matrix"
        os.system(cmd2)

    if (min(list(logomatI.index)) == 0):
        logomatI.index = logomatI.index + 1
     
    if (pos is not None):
        positions = [int(i) for i in pos]
        for posTmp in positions:
            logomatI.loc[posTmp]=0
            
    return logomatI
    
def main():
    args = parse_arguments()
    refArgs = vars(args)
    
    myFiles = refArgs.get("f") 
    myTitles = refArgs.get("titles") 
    myLogoType = refArgs.get("type") ## logo or ice logo
    optionInput = refArgs.get("optionInput") ##input file type - peptide list or a pssm matrix 
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

    
    if (len(myTitles) != len(myFiles)):
        for i in range(0, len(myFiles)):
            myTitles.append("")
    
    nrows = len(myFiles)
    ncols = 1
    height = 2.5
    fig = plt.figure(figsize=(10,height * nrows))
    for i in range(0,len(myFiles)):
        myFile = myFiles[i]
        myTitle = myTitles[i]
        
        logomatI = getLogomatI(optionInput, pos, myLogoType, myFile)
        axis = fig.add_subplot(nrows,ncols,i+1)
        showlogo(logomatI,myTitle,myLogoType, axis)
    
    
    fig.tight_layout()
    plt.savefig("myLogo.pdf")
    plt.close('all')
    
    

if __name__ == '__main__':
    main()



