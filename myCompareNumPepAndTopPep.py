#!/bin/python
import argparse
#import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f',  nargs='*', 
                                    help='The input should be lists of peptides ')
   
    my_parser.add_argument('-title', nargs='?',default='', help='title of the plot')
    my_parser.add_argument('-rounds', nargs='*',help='rounds of the plots')
    my_parser.add_argument('-o', nargs='?',default='TopNumber', help='Top * will be selected')
    my_parser.add_argument('-largest', nargs='?',default='10',help='**largest values of the data')
    
    return my_parser.parse_args()


def plot_Compare(df, myTitle, yLabel):
    fig,ax = plt.subplots()
    ax.plot(df.index, df['uniqCount'], color="red", marker="o",label='Unique peptides')
    #ax.set_xlabel("year",fontsize=14)
    ax.set_ylabel("Unique peptides",color="red",fontsize=14)
    # twin object for two different y-axis on the sample plot
    ax2=ax.twinx()
    # make a plot with different y-axis using second axis object
    ax2.plot(df.index, df['percTop'],color="blue",marker="o",label=yLabel)
    ax2.set_ylabel("% Total reads",color="blue",fontsize=14)

    ax.legend(loc='center left', frameon=False, bbox_to_anchor=(-0.05, -0.15))
    ax2.legend(loc='center right', frameon=False, bbox_to_anchor=(1., -0.15))

    ax.set_title(myTitle, pad=10)

    #fig.suptitle(myTitle, fontsize=20)

    fig.tight_layout()
    fig.savefig("CompareNumPepAndTopPep.pdf")
    plt.close("all")

def main():
    
    args = parse_arguments()
    refArgs = vars(args)
    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFiles=refArgs.get("f") ##rounds file
    myRounds = refArgs.get("rounds") ## rounds of the plots
    myTitle = refArgs.get("title") ## title of the plot
    option = refArgs.get("o") ## Top * will be selected
    largest = int(refArgs.get("largest") ) ## **largest values of the data
    
    yLabel = ""
    
    if (len(myRounds) != len(myFiles) or "" in myRounds):
        sys.exit("Please provide a round name for each input file, e.g. R0, R1 or R2.")
         
    refKey = {}
    refKey["Rounds"] = ["uniqCount", "percTop"]
    for i in range(0,len(myFiles)):

        cmd0 = ''' grep -v "#" ''' + myFiles[i] + " > myTmp0"
        os.system(cmd0)
        df_Round = pd.read_csv("myTmp0",header=0,index_col=0,sep="\t")
        uniqCount = len(df_Round.index)
        myRound = myRounds[i]    
        totalReads = df_Round['Total_Count'].sum()

        if (option == "TopNumber"):
            df_RoundTop =(df_Round.nlargest(largest, ['Frequency']))
            yLabel = 'Top ' + str(largest) + ' peptides (% of total reads)'
        else:
            df_RoundTop =(df_Round.nlargest(int(largest * uniqCount /100.0), ['Frequency']))
            yLabel = 'Top ' + str(largest) + '% peptides (% of total reads)'
        
        #print (df_RoundTop.shape)
        totalReadsTop = df_RoundTop['Total_Count'].sum()

        percTop = totalReadsTop / totalReads * 100    
        refKey[myRound] = [uniqCount, percTop]

        cmd1 = "rm myTmp0"
        os.system(cmd1)

    
    df = pd.DataFrame.from_dict(refKey, orient="index")
    new_header = df.iloc[0] 
    df = df[1:] 
    df.columns = new_header

    df = df.sort_index()
    plot_Compare(df, myTitle, yLabel)

    
if __name__ == '__main__':
    main()

