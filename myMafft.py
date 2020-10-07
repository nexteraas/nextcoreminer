#!/bin/env python3
import argparse
import os
import sys
import pandas as pd
import numpy as np
from math import ceil


def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='?',  help='Peptide file ')
    my_parser.add_argument('-op', nargs='?',  default="1.53", help='Gap opening penalty')
    my_parser.add_argument('-ep', nargs='?',  default="0", help='Gap extension penalty')
    my_parser.add_argument('-o', nargs='?', default="fasta",help='fasta or clustal output')
    

    return my_parser.parse_args()


def main():
    args = parse_arguments()
    refArgs = vars(args)

    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFile=refArgs.get("f") ##peptide file
    op = refArgs.get("op").strip() ##Gap opening penalty
    ep = refArgs.get("ep").strip() ##gap extension penalty
    option = refArgs.get("o") ##fasta or clustal output
    clustal = ""
    if (option == "clustal"):
        clustal = " --clustalout"
    
    
    cmd0 = "cut -f 1 " + myFile + ''' |awk '$0 ~ /^[A-Z]+$/' |gawk '{print ">"$0"\\n"$0}' > myTmp0.fasta'''
    os.system(cmd0)
    cmd1 = "mafft --quiet --thread 8  --reorder --retree 2 --op " + op + " --ep " + ep  + clustal + " --auto myTmp0.fasta > myMafft.txt"
    os.system(cmd1)

    cmd2 = "rm myTmp0.fasta"
    os.system(cmd2)

if __name__ == '__main__':
    main()



