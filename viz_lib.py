#!/usr/bin/python3
import pandas as pd
import numpy as np
import re

import matplotlib
# Turn interactive plotting off, set new render
matplotlib.use('Agg')
import matplotlib.pyplot as plt




def form_dict_to_df(dict_in):
	try:
		data_df = pd.DataFrame(dict_in.items(), columns = ['ref','value'])
		data_df['value'] = data_df['value'].astype(int)
		return data_df
	except TypeError:
		print("the object passed was not a dict")

def check_barcode(barcode):
	if re.search(r'^[ATCG]*$', barcode):
		return 1
	else:
		return 0

def get_max_value(*args):
	max = 0
	if len(args) > 1:
		for data in args:    
			tmp_max = get_max_value(data['value'])
			
			if max < tmp_max:
				max = tmp_max
		return max
	else:
		return args[0].max()

def data_overview_viz(name, unique_sumstats, retention_sumstats, reading_frame, barcodes):

	unique_sumstats_df = form_dict_to_df(unique_sumstats)
	retention_sumstats_df = form_dict_to_df(retention_sumstats)
	reading_frame_df = form_dict_to_df(reading_frame)
	barcodes_df = form_dict_to_df(barcodes)


	barcodes_df['selection'] = barcodes_df['ref'].apply(check_barcode)

	barcodes_df['rank'] = barcodes_df['value'].rank(ascending = True)
	barcodes_df = barcodes_df.sort_values(by = 'rank')

	x_limit_value = get_max_value(unique_sumstats_df, retention_sumstats_df, reading_frame_df, barcodes_df)

	x_limit_value= 1.05 * x_limit_value
	reading_frame_df = reading_frame_df.sort_values(by = 'ref')



	fig = plt.figure(figsize = (12,12))
	## Unique peptide data oveview plot
	ax = fig.add_subplot(2,2,1)
	ax.grid(False)
	ax.set_ylim(0,x_limit_value)
	ax.set_ylabel('Count')
	ax.set_title('Count of unique sequences at each step', fontsize = 20)

	x_span = len(unique_sumstats_df['ref'])
	x_post = np.arange(0,x_span,1)
	w = 0.98
	plt.xticks(x_post, unique_sumstats_df['ref'], rotation=10, horizontalalignment='right')
	plt.bar(x_post, unique_sumstats_df['value'], width = w, align = 'center', ec = "#82b8d0",fc = "#4f7282")

	## total count sequnce, and data retention plot
	ax = fig.add_subplot(2,2,2)
	ax.grid(False)
	ax.set_ylim(0,x_limit_value)
	ax.set_ylabel('Count')
	ax.set_title('Sequence count at each step', fontsize = 20)

	x_span = len(retention_sumstats_df['ref'])
	x_post = np.arange(0,x_span,1)
	w = 0.98
	plt.xticks(x_post, retention_sumstats_df['ref'], rotation=10, horizontalalignment='right')
	plt.bar(x_post, retention_sumstats_df['value'], width = w, align = 'center', ec = "#82b8d0",fc = "#4f7282")

	## Readning frame abundance plot
	ax = fig.add_subplot(2,2,3)

	max_val = int(reading_frame_df['ref'].max())
	axis_range = list(range(3,33,3))
	axis_range.extend(range(33, max_val+3 ,3))

	ax.grid(False)
	ax.set_yticks(axis_range)
	ax.set_xlim(0,x_limit_value)
	ax.set_xlabel('Count')
	ax.set_title('Kmer abundance', fontsize = 20)
	plt.barh(reading_frame_df['ref'], reading_frame_df['value'], align = 'center', ec = "#82b8d0",fc = "#4f7282")

	## Barode abundance overview plot
	ax = fig.add_subplot(2,2,4)
	ax.grid(False)
	ax.set_xlim(0,x_limit_value)
	ax.set_xlabel('Count')
	ax.set_title('Barcode abundance', fontsize = 20)
	plt.barh(barcodes_df[barcodes_df['selection'] == 1]['ref'], barcodes_df[barcodes_df['selection'] == 1]['value'], align = 'center', ec = "#82b8d0",fc = "#4f7282")

	plt.savefig(name + '.sum_stats.png')
	plt.close(fig)
