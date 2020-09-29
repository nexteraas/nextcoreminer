#!/bin/env python3
import argparse
import os
import sys
import pandas as pd


def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='?',  help='The matrix file from myVenn.py ')

    my_parser.add_argument('-g',  nargs='*',
                                    help='The input should be the groups of interest')


    return my_parser.parse_args()

    
#### given the groups, retrive the peptides: 
def readFile (myFile, myGroups):
    file0 = open(myFile)
    out =  open("myPepFromVenn.txt","w")
    
    while True:
        line = file0.readline().rstrip()
        if (len(line)==0):
            break
        
        arr = line.split("\t")
        if (arr[1] == "group" or arr[1] in myGroups):
            out.write(line + "\n")
        
    file0.close()
    out.close()

def main():
    args = parse_arguments()
    refArgs = vars(args)

    if ( "f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFile = refArgs.get("f") ##The matrix file from myVenn.py
    myGroups = refArgs.get("g") ##groups of interest

    readFile(myFile,myGroups)

if __name__ == "__main__":
    main()
