#!/usr/bin/env python

import sys
import gzip
import numpy
import pyBigWig

from pathlib import Path

def main():
	if len(sys.argv) < 6:
		print ("Usage : " + sys.argv[0] + " [bin_of_interest_bedgz_filename] [bin_size] [tile_size] [signal_bigwig_filename] [matrix_npy_filename]")
		sys.exit()
	else:
		bin_of_interest_bedgz_filename = sys.argv[1]
		bin_size = int(sys.argv[2])
		tile_size = int(sys.argv[3])
		signal_bigwig_filename = sys.argv[4]
		matrix_npy_filename = sys.argv[5]
		build_tiled_epigenomic_feature(bin_of_interest_bedgz_filename,bin_size,tile_size,signal_bigwig_filename,matrix_npy_filename)


def build_tiled_epigenomic_feature(bin_of_interest_bedgz_filename,bin_size,tile_size,signal_bigwig_filename,matrix_npy_filename):
	bin_of_interest_list = load_bin_of_interest(bin_of_interest_bedgz_filename)
	print(tile_size)
	tiled_epigenome_signal_matrix = extract_tiled_epigenome_matrix(bin_of_interest_list,bin_size,tile_size,signal_bigwig_filename)
	Path.mkdir(Path(matrix_npy_filename).parent, parents=True, exist_ok=True)
	numpy.save(matrix_npy_filename,tiled_epigenome_signal_matrix)


def load_bin_of_interest(bin_of_interest_bedgz_filename):
	bin_of_interest_list = []
	with gzip.open(bin_of_interest_bedgz_filename,"rt") as bin_of_interest_bedgz_file:
		for rawline in bin_of_interest_bedgz_file:
			fields = rawline.strip().split()
			bin_of_interest = (fields[0],int(fields[1]),int(fields[2]))
			bin_of_interest_list.append(bin_of_interest)

	return bin_of_interest_list


def extract_tiled_epigenome_matrix(bin_of_interest_list,bin_size,tile_size,signal_bigwig_filename):
	tile_count = int(bin_size/tile_size)
	tiled_epigenome_signal_matrix = numpy.zeros((len(bin_of_interest_list),tile_count),dtype=numpy.uint16)
	with pyBigWig.open(signal_bigwig_filename) as signal_bigwig_file:
		for bin_index,bin_of_interest in enumerate(bin_of_interest_list):
			try:
				epigenome_signal_vector = list(map(integize,numpy.nan_to_num(numpy.array(list(map(remove_none,signal_bigwig_file.stats(*bin_of_interest,type="sum",nBins=tile_count))))*1/tile_count)))
				tiled_epigenome_signal_matrix[bin_index,:] = epigenome_signal_vector
			except:
				print (str(bin_index) + str(bin_of_interest))

	return tiled_epigenome_signal_matrix

def remove_none(tiled_epigenome_signal):
	return 0 if tiled_epigenome_signal is None else tiled_epigenome_signal 


def integize(base_signal):
	# adjusted_signal = base_signal*100 if base_signal*100 <= 65535 else  65535 #  original 
	adjusted_signal = base_signal*1000 if base_signal*1000 <= 65535 else  65535
	return numpy.uint16(adjusted_signal)
	


if __name__ == "__main__":
	main()



