#!/bin/python3
import pandas as pd
def getDynamicOptions(fileIn):
    options = ["te","test"]
    df = pd.read_csv(fileIn,index_col=0, header=0,sep="\t",dtype=str)      
    return (options)
    #return (df.group.unique())
#getDynamicOptions ("dataset_2934.dat")


