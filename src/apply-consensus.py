#! /usr/bin/env python

# Copyright Ivan Sovic, 2015. www.sovic.org
#
# Applies the variants found by consensus.py script to a given reference.

import os;
import sys;
import operator;
import subprocess;

if __name__ == "__main__":
	if (len(sys.argv) < 5):
		sys.stderr.write('Applies the variants found by consensus.py script to a given reference.\n');
		sys.stderr.write('\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <reference_file_path> <variant_file> <{sb}am_file_> [position]\n' % sys.argv[0]);
		sys.stderr.write('\t(If <collective_output_file> is equal to "-", no files will be written to disk.)\n');
		sys.stderr.write('\tPosition parameter is a string specifying "chromosome:start-end"\n\n');
		exit(1);
	
	reference_file = sys.argv[1];
	coverage_threshold = int(sys.argv[2]);
	output_prefix = sys.argv[3];
	sam_file = sys.argv[4];
	bed_position = '';
	if (len(sys.argv) > 5):
		bed_position = sys.argv[5];
	# sys.stderr.write('bed_position: "%s"\n\n' % bed_position);
	
	processes = [];

	if (output_prefix == '-'):
		output_prefix = os.path.splitext(sam_file)[0];
	main(sam_file, reference_file, coverage_threshold, output_prefix, 0, bed_position);

	if (output_prefix != '-'):
		CollectSummaries([sam_file], output_prefix + '.variant.sum');
