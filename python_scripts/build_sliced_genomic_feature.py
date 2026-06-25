#!/usr/bin/env python

import sys
import gzip
import itertools 
import numpy 
from pathlib import Path

def main():
	if len(sys.argv) < 5:
		print ("Usage : " + sys.argv[0] + " [anchor_label_txtgz_filename] [chrom] [seq_numpy_filename] [sliced_seq_array_filename]")
		sys.exit()
	else:
		anchor_label_txtgz_filename = sys.argv[1]
		chrom = sys.argv[2]
		seq_numpy_filename = sys.argv[3]
		sliced_seq_array_filename = sys.argv[4]
		build_sliced_genomic_feature(anchor_label_txtgz_filename , chrom  , seq_numpy_filename , sliced_seq_array_filename)


def build_sliced_genomic_feature(anchor_label_txtgz_filename , chrom  , seq_numpy_filename , sliced_seq_array_filename):
	anchor_range = extract_chrom_anchor_range(anchor_label_txtgz_filename,chrom)
	sliced_seq_array = load_sliced_seq_array(anchor_range,seq_numpy_filename)
	Path.mkdir(Path(sliced_seq_array_filename).parent, parents=True, exist_ok=True)
	numpy.save(sliced_seq_array_filename,sliced_seq_array)


def load_sliced_seq_array(anchor_range,seq_numpy_filename):
	min_index,max_index = anchor_range
	seq_array = numpy.load(seq_numpy_filename,mmap_mode="r")
	sliced_seq_array = numpy.unpackbits(seq_array[min_index:max_index+1,],axis=1)
	return sliced_seq_array



def load_epi_array(epi_numpy_filename ,anchor_range, cap_crit = 0.95):
	min_index,max_index = anchor_range
	epi_array = numpy.load(epi_numpy_filename,mmap_mode="r")
	cap_val = int(numpy.quantile(epi_array,cap_crit))
	epi_array = epi_array[min_index:max_index+1,]
	epi_array = numpy.clip(epi_array,None,cap_val) * 0.001
	epi_array = epi_array.astype(numpy.float32)
	return epi_array


def extract_chrom_anchor_range(anchor_label_txtgz_filename,chrom):
	index_list = []
	with gzip.open(anchor_label_txtgz_filename,"rt") as anchor_label_txtgz_file:
		for rawline in itertools.islice(anchor_label_txtgz_file,1,None):
			fields = rawline.strip().split()
			if fields[1].split(":")[0] == chrom:
				index = int(fields[0]) - 1
				index_list.append(index)

	anchor_range = (min(index_list),max(index_list))

	return anchor_range


if __name__ == "__main__":
	main()


