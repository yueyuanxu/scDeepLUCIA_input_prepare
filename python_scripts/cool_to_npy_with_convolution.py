#!/usr/bin/env python
import sys
import numpy
import cooler

from scipy.signal import convolve2d
from pathlib import Path

def main():
	if len(sys.argv) < 5:
		print ("Usage : " + sys.argv[0] + " [cool_filename] [chrom] [window_size] [convolved_contact_matrix_npy_filename]")
		sys.exit()
	else:
		cool_filename = sys.argv[1]
		chrom = sys.argv[2]
		window_size = int(sys.argv[3])
		convolved_contact_matrix_npy_filename = sys.argv[4]
		cool_to_convolved_npy(cool_filename , chrom , window_size , convolved_contact_matrix_npy_filename)


def cool_to_convolved_npy(cool_filename , chrom , window_size , convolved_contact_matrix_npy_filename):
	contact_matrix = cooler.Cooler(cool_filename).matrix(balance=False).fetch(chrom)
	kernel = numpy.ones((window_size,window_size)).astype(numpy.int32)
	convolved_matrix = convolve2d(contact_matrix, kernel, mode='same',boundary="symm")
	Path.mkdir(Path(convolved_contact_matrix_npy_filename).parent, parents=True, exist_ok=True)
	numpy.save(convolved_contact_matrix_npy_filename,convolved_matrix)


if __name__ == "__main__" :
	main()

