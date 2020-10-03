#!/bin/env python3
import argparse
import os
import sys
import pandas as pd
import re

def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='*',  help='The clusters of peptides ')

    my_parser.add_argument('-peps',  nargs='?',
                                    help='The input should be the peptides of interest.')


    return my_parser.parse_args()

    

def readFile (myFile):
    file0 = open(myFile)
    refTmpPeps= {}
    cluster=""
    while True:
        line = file0.readline().rstrip()
        if (len(line)==0):
            break
        arr = line.split("\t")
        refTmpPeps[arr[0]] = arr[1] ##not arr[1:]
        if ("Cluster" in arr[0]):
            cluster = arr[0].split("_")[0]
    file0.close()
    return (refTmpPeps,cluster)

def main():
    args = parse_arguments()
    refArgs = vars(args)

    if ( "f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFiles= refArgs.get("f") ##clusters of peptides
    myPep = refArgs.get("peps") ##peptides of interest
    delimiters = ["X",",", ";", " "]
    regexPattern = '|'.join(map(re.escape, delimiters))
    myPepTmp= re.split(regexPattern, myPep)
    myPeps = [x for x in myPepTmp if x]


    out = open("myExactPeps.txt","w")

    refExactPeps={}
    for tmpFile in myFiles:
        refPeps,cluster = readFile(tmpFile)

        for pep in myPeps:
            if (pep in refPeps):
                refExactPeps[pep] = cluster + "\t" + refPeps.get(pep) 

    out.write("Peptides" + "\t" + "Clusters" + "\t" + "group" + "\n")
    for keys in refExactPeps:
        out.write(keys + "\t" + refExactPeps.get(keys)  + "\n")
    out.close()

if __name__ == "__main__":
    main()
