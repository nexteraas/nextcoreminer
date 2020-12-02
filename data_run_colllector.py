#!/usr/bin/python3 
import argparse
import os
import re
import sys

import datetime as dt

import pandas as pd
import json

from pathlib import Path
import string
import random

import concurrent.futures

import data_preprocessign_methods as dpm
import file_handler as fh

def parse_arguments():
	""" All arugments must be there in oder for the program to run"""

	parser = argparse.ArgumentParser()
	
	parser.add_argument('-d', '--datafiles', default=sys.stdin, required=True, nargs='+',
						help='Path to the input files'
						'The output must be from the data_cleaing.py scirpt in order to work'
						)
	parser.add_argument('-p', '--processes', default=4, required=True,
						help='Numver of processes to be used when reading in the files'
						'Be aware of the max number of possible on the system'
						)
	parser.add_argument('-t', '--topfolder', required=True,
						help='Number of processes to be used when reading in the files'
						)
	return parser.parse_args()

def generate_hyperlink(web_folder):
 	"""Generate the clicable hyperlink to be displayed in galaxy"""
 	with open('hyperlink.html', 'w') as of:
 		of.write("<!DOCTYPE html>")
 		of.write("<html>")
 		of.write("<head>")
 		of.write("<title>Nextera Data table hyperlink</title>")
 		of.write("</head>")
 		of.write("<body>")
 		of.write("<a href=\"http://172.16.10.2:5000/#/"+web_folder+"/\" target=\"_blank\">Go to data location</a>")
 		of.write("</body>")
 		of.write("</html>")

def main():
	args = parse_arguments()
	if args:
		
		list_peptide_data = fh.data_list_handler_json(args.datafiles, int(args.processes))
		
		dict_peptidets_per_run, list_rounds, NNK_name = dpm.collect_peptide_dict_from_json(list_peptide_data)
		
		work_folder, web_folder= fh.generate_workfolder(NNK_name, args.topfolder)
		
		fh.data_for_json(dict_peptidets_per_run, list_rounds, work_folder)

		dpm.peptide_nice_fromat(dict_peptidets_per_run, list_rounds)
		generate_hyperlink(web_folder)		

if __name__ == '__main__':
	main()

