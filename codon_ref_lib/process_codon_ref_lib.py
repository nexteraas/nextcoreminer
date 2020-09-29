#!/usr/bin/python3
import argparse
import os
import re
import sys
import pandas as pd
import numpy as np
import json

def parse_arguments():
	""" All arugments must be there in oder for the program to run"""

	parser = argparse.ArgumentParser()
	
	parser.add_argument('-l', '--libfile', default=sys.stdin, required=True,
						help='Path to the libray file containing the condon reference libraries'
						)

	return parser.parse_args()


def process_lib_file(lib_file_in):
	codon_dict = None
	codon_table_name = None
	with open(lib_file_in, 'r') as f:
		
		for line in f:

			_line = line.rstrip('\r\n')

			if re.search(r'^>', _line) is not None:
				
				if codon_dict is not None:
					codon_json = json.dumps(codon_dict)
					of = open(codon_table_name+".json","w")
					of.write(codon_json)
					of.close()

				
				codon_table_name = _line[1:]
				codon_table_name = codon_table_name.lower()
				print(codon_table_name)
				codon_dict = {}
			else:
				codon, aa = _line.split('\t')
				codon_dict[codon] = aa


def main():
	args = parse_arguments()

	if args:
		process_lib_file(args.libfile)




if __name__ == "__main__":
	main()