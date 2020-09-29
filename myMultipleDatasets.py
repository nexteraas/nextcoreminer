import pandas as pd
import sys
import argparse
import os

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f',  nargs='*', 
                                    help='The input should be lists of peptides with '
                                    'specificity and frequency ')
    return my_parser.parse_args()

def mergeDF(myFiles):
    df_merged = pd.DataFrame()
    refCol = {}
    
    for files_i in myFiles:
        
        cmd0 =  ''' gawk -F "\t" 'length($1)==0 ||$1=="AA_Peptide" || $1 ~ /^[A-Z]+$/' ''' + files_i  + "> myTmp0"  
        os.system(cmd0)
        df0 = pd.read_csv("myTmp0",header=0,index_col=0,sep="\t")
        
        columns = list(df0.columns)
        for col in columns:
            refCol[col] = 1
        df_merged = pd.concat([df_merged,df0],axis=0,sort=True)
        cmd1 = "rm myTmp0"
        os.system(cmd1)
        
    sortedCols = list(refCol)
    df_merged = df_merged[sortedCols]
    
    df_merged = df_merged[~df_merged.index.duplicated(keep='first')]
   
    df_merged.fillna("NA").to_csv("merged_matrix.txt", sep="\t", index = True, header=True)
    
def main():  
    args = parse_arguments()
    refArgs = vars(args)
    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")
     
    myFiles=refArgs.get("f") ##rounds file
    mergeDF(myFiles)

if __name__ == '__main__':
    main()


