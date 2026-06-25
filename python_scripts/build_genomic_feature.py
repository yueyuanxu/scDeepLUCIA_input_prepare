#!/usr/bin/env python

import sys
import gzip
import numpy
import pysam

from pathlib import Path

def main():
	if len(sys.argv) < 5:
		print ("Usage : " + sys.argv[0] + " [bin_of_interest_bedgz_filename] [bin_size] [genome_fasta_filename] [matrix_npy_filename]")
		sys.exit()
	else:
		bin_of_interest_bedgz_filename = sys.argv[1]
		bin_size = int(sys.argv[2])
		genome_fasta_filename = sys.argv[3]
		matrix_npy_filename = sys.argv[4]
		build_genomic_feature(bin_of_interest_bedgz_filename,bin_size,genome_fasta_filename,matrix_npy_filename)


def build_genomic_feature(bin_of_interest_bedgz_filename,bin_size,genome_fasta_filename,matrix_npy_filename):
	bin_of_interest_list = load_bin_of_interest(bin_of_interest_bedgz_filename)
	one_hot_encoding_matrix = extract_one_hot_encoding_matrix(bin_of_interest_list,bin_size,genome_fasta_filename)
	one_hot_encoding_matrix = numpy.packbits(one_hot_encoding_matrix,1)
	Path.mkdir(Path(matrix_npy_filename).parent, parents=True, exist_ok=True)
	numpy.save(matrix_npy_filename,one_hot_encoding_matrix)


def load_bin_of_interest(bin_of_interest_bedgz_filename):
	bin_of_interest_list = []
	with gzip.open(bin_of_interest_bedgz_filename,"rt") as bin_of_interest_bedgz_file:
		for rawline in bin_of_interest_bedgz_file:
			fields = rawline.strip().split()
			bin_of_interest = (fields[0],int(fields[1]),int(fields[2]))
			bin_of_interest_list.append(bin_of_interest)

	return bin_of_interest_list


def extract_one_hot_encoding_matrix(bin_of_interest_list,bin_size,genome_fasta_filename):
	genome_fasta = pysam.Fastafile(genome_fasta_filename)
	one_hot_encoding_matrix = numpy.zeros((len(bin_of_interest_list),bin_size,4),dtype=numpy.uint8)
	for bin_index,bin_of_interest in enumerate(bin_of_interest_list):
		bin_seq = genome_fasta.fetch(*bin_of_interest).upper()
		for base_index,base in enumerate(bin_seq):
			one_hot_encoding_matrix[bin_index,base_index,:] = encode(base)

	return one_hot_encoding_matrix


def encode(base):
	base_to_encode = {
		"A": numpy.asarray([1,0,0,0],dtype=numpy.uint8),
		"C": numpy.asarray([0,1,0,0],dtype=numpy.uint8),
		"G": numpy.asarray([0,0,1,0],dtype=numpy.uint8),
		"T": numpy.asarray([0,0,0,1],dtype=numpy.uint8),
		"N": numpy.asarray([1,1,1,1],dtype=numpy.uint8)}
	return base_to_encode[base]



if __name__ == "__main__":
	main()
