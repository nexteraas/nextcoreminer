#!/bin/env python3
import argparse
import os
import sys
import pandas as pd
import re

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='?',  help='The matrix file from myVenn.py ')

    my_parser.add_argument('-g',  nargs='?',
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

        if (arr[1] == "group" ):
            out.write(line + "\n")
        else:
            tmp = arr[1].split("+")
            tmp.sort()
            arr[1] = "+".join(tmp)
            
            if (arr[1] in myGroups):
                out.write(line + "\n")

    file0.close()
    out.close()

def main():
    args = parse_arguments()
    refArgs = vars(args)

    if ( "f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFile = refArgs.get("f") ##The matrix file from myVenn.py
    myGroup = refArgs.get("g") ##groups of interest
    delimiters = ["X",",", ";", " "]
    regexPattern = '|'.join(map(re.escape, delimiters))
    myGroupTmp= re.split(regexPattern, myGroup)
    myGroups = [x for x in myGroupTmp if x]
    myGroups_sorted = []

    for myGroup in myGroups:
        tmp = myGroup.split("+")
        tmp.sort()
        myGroups_sorted.append("+".join(tmp))


    readFile(myFile,myGroups_sorted)

if __name__ == "__main__":
    main()
