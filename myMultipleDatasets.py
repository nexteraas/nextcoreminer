import pandas as pd
import sys
import argparse
import os

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-fInfoOption',  nargs='?',default="no", help="Do you have a file with meta-data of peptides")
    
    my_parser.add_argument('-fInfo',  nargs='?', 
                                    help='The input should be lists of peptides with counts frequency ')
    
    my_parser.add_argument('-f',  nargs='*', 
                                    help='The input should be lists of peptides ')
    return my_parser.parse_args()

def mergeDF(myFiles, myFileInfo,fInfoOption):
    df_merged = pd.DataFrame()
    refCol = {}
    
    for files_i in myFiles:
        
        cmd0 =  ''' gawk -F "\t" '{if (NR==1 && $1~ /^[A-Z]+$/) print "AA_Peptide\\n"$0;else print $0 }' ''' + files_i  + "> myTmp0"  
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
     
    
    if (fInfoOption == "yes"):
        dfInfo = pd.read_csv(myFileInfo,header=0,index_col=0,sep="\t")

        dfInfo = dfInfo[dfInfo.index.isin (df_merged.index)]

        columns = list(dfInfo.columns)
        for col in columns:
                refCol[col] = 1


        sortedCols = list(refCol)

        df_merged = pd.concat([dfInfo,df_merged],axis=0,sort=True)

    df_merged = df_merged[~df_merged.index.duplicated(keep='first')]
    df_merged = df_merged[sortedCols]

    df_merged = df_merged.rename_axis('Peptides')
    df_merged.fillna("NA").to_csv("merged_matrix.txt", sep="\t", index = True, header=True)
    
    df_merged2 =pd.DataFrame(df_merged.to_records())
    df_merged2.fillna("NA").to_json("merged_matrix.json")
    
   

def main():  
    args = parse_arguments()
    refArgs = vars(args)
    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")
     
    myFiles = refArgs.get("f") 
    myFileInfo = refArgs.get("fInfo") 
    myFileInfoOption = refArgs.get("fInfoOption")
    mergeDF(myFiles,myFileInfo,myFileInfoOption)

if __name__ == '__main__':
    main()



