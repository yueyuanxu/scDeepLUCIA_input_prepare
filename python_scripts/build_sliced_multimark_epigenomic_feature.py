#!/usr/bin/env python

import sys
import gzip
import itertools 

import numpy

from pathlib import Path


def main():
	if len(sys.argv) < 6:
		print ("Usage : " + sys.argv[0] + " [bin_of_interest_bedgz_filename] [chrom] [epi_array_dirname] [mark_list] [multimark_epi_array_filename]")
		sys.exit()
	else:
		bin_of_interest_bedgz_filename = sys.argv[1]
		chrom = sys.argv[2]
		epi_array_dirname = sys.argv[3]
		mark_list = sys.argv[4].split(",")
		multimark_epi_array_filename = sys.argv[5]
		build_sliced_multimark_epigenomic_feature(bin_of_interest_bedgz_filename , chrom  , epi_array_dirname , mark_list, multimark_epi_array_filename)


def build_sliced_multimark_epigenomic_feature(bin_of_interest_bedgz_filename , chrom  , epi_array_dirname , mark_list, multimark_epi_array_filename):
	anchor_range = extract_chrom_anchor_range(bin_of_interest_bedgz_filename,chrom)
	mark_to_epi_array_filename = load_epi_array_dir(epi_array_dirname)
	multimark_epi_array = load_multimark_epi_array(mark_list,mark_to_epi_array_filename,anchor_range)
	Path.mkdir(Path(multimark_epi_array_filename).parent, parents=True, exist_ok=True)
	numpy.save(multimark_epi_array_filename,multimark_epi_array)


def load_multimark_epi_array(mark_list,mark_to_epi_array_filename,anchor_range):
	epi_array_list = []
	for mark in mark_list:
		epi_numpy_filename = mark_to_epi_array_filename[mark]
		epi_array = load_epi_array(epi_numpy_filename,anchor_range)
		epi_array_list.append(epi_array)
	multimark_epi_array = numpy.stack(epi_array_list,axis=2)
	return multimark_epi_array


def load_epi_array(epi_numpy_filename ,anchor_range, cap_crit = 0.95):
	min_index,max_index = anchor_range
	epi_array = numpy.load(epi_numpy_filename,mmap_mode="r")
	cap_val = int(numpy.quantile(epi_array,cap_crit))
	epi_array = epi_array[min_index:max_index+1,]
	epi_array = numpy.clip(epi_array,None,cap_val) * 0.001
	epi_array = epi_array.astype(numpy.float32)
	return epi_array


def extract_chrom_anchor_range(bin_of_interest_bedgz_filename,chrom):
	index_list = []
	with gzip.open(bin_of_interest_bedgz_filename,"rt") as bin_of_interest_bedgz_file:
		for index,rawline in enumerate(bin_of_interest_bedgz_file):
			fields = rawline.strip().split()
			if fields[0] == chrom:
				#index = int(fields[3].split("_")[1])
				#min_index = index if index < min_index else min_index 
				#max_index = index if max_index < index else max_index
				index_list.append(index)

	anchor_range = (min(index_list),max(index_list))

	return anchor_range


def load_epi_array_dir(epi_array_dirname):
	mark_to_epi_array_filename ={}
	epi_array_dir = Path(epi_array_dirname)
	for epi_array_filename  in epi_array_dir.glob("*.npy"):
		mark = epi_array_filename.name.split(".")[0]
		mark_to_epi_array_filename[mark] = epi_array_filename

	return mark_to_epi_array_filename


if __name__ == "__main__":
	main()


