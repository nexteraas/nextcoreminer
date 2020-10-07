#!/bin/env python3
import argparse
import os
import sys
import pandas as pd
import numpy as np
from math import ceil


def parse_arguments():
    my_parser = argparse.ArgumentParser(allow_abbrev=False)

    my_parser.add_argument('-f', nargs='?',  help='The mafft alignment output/psiblast checkpoint file ')

    my_parser.add_argument('-db',  nargs='?', default="uniprot_sprot",
                                    help='Blast database ')
    my_parser.add_argument('-evalue', nargs='?', default="10",help='E-value cut off')
    my_parser.add_argument('-op', nargs='?',  default="32767", help='Gap opening penalty')
    my_parser.add_argument('-ep', nargs='?',  default="32767", help='Gap extension penalty')
    my_parser.add_argument('-num', nargs='?',  default="1", help='number of iterations')
    my_parser.add_argument('-o', nargs='?',  default="msa", help='input file option')

    return my_parser.parse_args()




def main():
    args = parse_arguments()
    refArgs = vars(args)

    if ("f" not in refArgs):
        sys.exit("Please provide a valid file.")

    myFile=refArgs.get("f") ##clusters of peptides
    mydb = refArgs.get("db") ##blast db
    evalue = refArgs.get("evalue") ## evalue threshold cutoff
    op = refArgs.get("op") ##Gap opening penalty
    ep = refArgs.get("ep") ##gap extension penalty
    num = refArgs.get("num") ##number of iterations
    option = refArgs.get("o") ##input file type - msa or psiblast pssm checkpoint

    blastdb = "/opt/nextera/galaxy/data/blastdb/" + mydb
    
    if (option == "msa"):
        cmd1 = "psiblast -in_msa  " +  myFile + " -db " + blastdb + " -num_threads 10  -comp_based_stats 0 -gapopen " + op + " -gapextend " + ep  + " -evalue " + evalue + " -num_iterations " + num + " -out myPsiBlast.out.txt -out_pssm myPsiblast.pssm.checkpoint.txt   -out_ascii_pssm myPsiblast.pssm.matrix.txt  -ignore_msa_master "
        os.system(cmd1)
    
    else:
        cmd2 = "psiblast -in_pssm  " + myFile + " -db " + blastdb + " -num_threads 10  -comp_based_stats 0 -gapopen " + op  + " -gapextend " + ep + " -evalue " + evalue + " -num_iterations " + num + " -out myPsiBlast.out.txt -out_pssm myPsiblast.pssm.checkpoint.txt -out_ascii_pssm myPsiblast.pssm.matrix.txt "
        os.system(cmd2)

if __name__ == '__main__':
    main()



