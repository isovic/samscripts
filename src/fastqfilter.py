#! /usr/bin/python

# Copyright Ivan Sovic, 2015. www.sovic.org
#
# A script that takes an input FASTA/FASTQ file and performs various methods
# of filtering on it. Examples include: selecting sequences that contain certain
# strings in their header, and others.

import re;
import os;
import sys;
import fastqparser;
import numpy as np;



def filter_seqs_by_header(input_fastq_path, header_patterns_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	try:
		fp_filter = open(header_patterns_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % header_patterns_path);
		exit(0);

	filter_headers = fp_filter.readlines();
	filter_headers = [line.strip() for line in filter_headers if (len(line.strip()) > 0)];
	fp_filter.close();

	num_matches = 0;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			break;
		
		for filter_header in filter_headers:
			if ((filter_header[-1] == ';' and filter_header[0:-1].lower() == header.lower()) or
				(filter_header[-1] != ';' and filter_header.lower() in header.lower())):
				num_matches += 1;
				sys.stderr.write('\rFound %d seqs, last: "%s".' % (num_matches, header));
				fp_out.write('\n'.join(read) + '\n');
		
	sys.stderr.write('\n');
	fp_in.close();

def filter_duplicate_ncbi_id(input_fastq_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	num_matches = 0;
	id_hash = {};

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		header_id = header.split()[0].lower();

		if (header_id in id_hash):
			continue;

		num_matches += 1;
		sys.stderr.write('\rFound %d seqs, last: "%s".' % (num_matches, header));
		fp_out.write('\n'.join(read) + '\n');
		id_hash[header_id] = True;

	sys.stderr.write('\n');
	fp_in.close();

def join_fastq_lines(input_fastq_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	num_matches = 0;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();

def filter_seqs_by_length(input_fastq_path, min_length, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	num_matches = 0;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			break;
		
		if (len(read[1]) >= min_length):
			num_matches += 1;
			sys.stderr.write('\rFound %d seqs, last: "%s".' % (num_matches, header));
			fp_out.write('\n'.join(read) + '\n');
		# fp_out.write('\n');
		# fp_out.write('(1) ' + read[0] + '\n');
		# fp_out.write('(2) ' + read[1] + '\n');

			# exit(1);
		
	sys.stderr.write('\n');
	fp_in.close();

def base_quality_stats(input_fastq_path):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	num_reads = 0;
	means_1d = [];
	means_2d = [];

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		if (len(read) == 0):
			break;

		if (len(read) < 4):
			sys.stderr.write('Given file is not a FASTQ file! Exiting.\n');
			exit(1);

		quals = read[3];
		phreds = [];
		for char in quals:
			phreds.append(ord(char) - 33);
		if (('twodir' in header.lower()) or ('-2d' in header.lower())):
			means_2d.append(np.mean(phreds));
			# sys.stdout.write('[2d] ');
		else:
			means_1d.append(np.mean(phreds));
		# 	sys.stdout.write('     ');
		# sys.stdout.write('avg = %5.2f\tstd = %5.2f\tmin = %2d\tmax = %2d\n' % (np.mean(phreds), np.std(phreds), min(phreds), max(phreds)));

		num_reads += 1;

		# if (len(read[1]) >= min_length):
		# 	num_matches += 1;
		# 	sys.stderr.write('\rFound %d seqs, last: "%s".' % (num_matches, header));
		# 	fp_out.write('\n'.join(read) + '\n');
		# fp_out.write('\n');
		# fp_out.write('(1) ' + read[0] + '\n');
		# fp_out.write('(2) ' + read[1] + '\n');

			# exit(1);
		
	sys.stdout.write('[1d] avg = %5.2f\tstd = %5.2f\tmin = %2d\tmax = %2d\n' % (np.mean(means_1d), np.std(means_1d), min(means_1d), max(means_1d)));
	sys.stdout.write('[2d] avg = %5.2f\tstd = %5.2f\tmin = %2d\tmax = %2d\n' % (np.mean(means_2d), np.std(means_2d), min(means_2d), max(means_2d)));

	sys.stderr.write('\n');
	fp_in.close();

def reverse_complement_seqs(input_fastq_path):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	num_reads = 0;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		if (len(read) == 0):
			break;

		if (len(read) == 2):
			sys.stdout.write('%s\n' % read[0]);
			sys.stdout.write('%s\n' % fastqparser.revcomp_seq(read[1]));
		elif (len(read) == 4):
			sys.stdout.write('%s\n' % read[0]);
			sys.stdout.write('%s\n' % fastqparser.revcomp_seq(read[1]));
			sys.stdout.write('%s\n' % read[2]);
			sys.stdout.write('%s\n' % read[3][::-1]);

		num_reads += 1;

	sys.stderr.write('\n');
	fp_in.close();

def hard_clip(input_fastq_path, num_bases):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	num_reads = 0;

	sys.stderr.write('Trimming %dbp from the beginning and ending of each sequence.\n' % num_bases);
	sys.stderr.write('Starting to process file "%s".\n' % input_fastq_path);
	while True:
		num_reads += 1;
		sys.stderr.write('\rProcessing seq %d...' % num_reads);

		[header, read] = fastqparser.get_single_read(fp_in);
		if (len(read) == 0):
			break;

		if (len(read[1]) <= (2*num_bases)):
			# sys.stderr.write('Skipping, len(read[1]) = %d, num_bases = %d\n' % (len(read[1]), num_bases) );
			continue;

		if (len(read) == 2):
			clipped_seq = read[1][num_bases:-num_bases];

			sys.stdout.write('%s_clipped\n' % read[0]);
			sys.stdout.write('%s\n' % clipped_seq);
		elif (len(read) == 4):
			clipped_seq = read[1][num_bases:-num_bases];
			clipped_qual = read[3][num_bases:-num_bases];

			sys.stdout.write('%s_clipped\n' % read[0]);
			sys.stdout.write('%s\n' % clipped_seq);
			sys.stdout.write('%s\n' % read[2]);
			sys.stdout.write('%s\n' % clipped_qual);

	sys.stderr.write('\n');
	fp_in.close();

def remove_special_chars_from_headers(input_fastq_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	num_matches = 0;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		# read[0] = read[0][0] + read[0][1:].replace();
		read[0] = read[0][0] + re.sub('[^0-9a-zA-Z]', '_', read[0][1:]); # re.sub("[|:", "_", read[0][1:]);
		fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();

def filter_for_marginalign(input_fastq_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);
	if (fp_out == None):
		if (input_fastq_path == out_fastq_path):
			sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
			exit(0);
		try:
			fp_out = open(out_fastq_path, 'w');
		except:
			sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % out_fastq_path);
			exit(0);

	num_matches = 0;
	header_hash = {};

	while True:
		[header, read] = get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		# read[0] = read[0][0] + read[0][1:].replace();
		if (len(read[1]) <= 50000):
			# read[0] = read[0][0] + re.sub('[^0-9a-zA-Z]', '_', read[0][1:]); # re.sub("[|:", "_", read[0][1:]);
			new_header = read[0][0] + re.sub('[^0-9a-zA-Z]', '_', read[0][1:]); # re.sub("[|:", "_", read[0][1:]);
			header_hash[new_header[1:]] = read[0][1:];
			read[0] = new_header;
			fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();
	return header_hash;

def wrap_fastq_file(input_fastq_path, wrap_length, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	current_read = 0;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		current_read += 1;

		read[1] = re.sub("(.{%d})"%(wrap_length), "\\1\n", read[1], 0, re.DOTALL);	### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
		if (len(read) == 4):
			read[3] = re.sub("(.{%d})"%(wrap_length), "\\1\n", read[3], 0, re.DOTALL);	### Wrap the qual line, because DALIGNER has a 9998bp line len limit.
		fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();

def convert_to_uppercase(input_fastq_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	current_read = 0;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		current_read += 1;

		read[1] = read[1].upper();
		fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();

def enumerate_headers(input_fastq_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);

	current_read = 0;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		
		if (len(read) == 0):
			break;

		current_read += 1;

		# read[1] = read[1].upper();
		read[0] = read[0][0] + str(current_read);
		fp_out.write('\n'.join(read) + '\n');

	sys.stderr.write('\n');
	fp_in.close();

### Modifies the FASTQ headers to replace all special chars with '_'.
def keep_gi_header(input_fastq_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);
	num_matches = 0;
	# header_hash = {};
	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		if (len(read) == 0):
			break;
		read[0] = read[0].split()[0];
		fp_out.write('\n'.join(read) + '\n');
	sys.stderr.write('\n');
	fp_in.close();
	fp_out.close();
	# return header_hash;

### Modifies the FASTQ headers to replace all special chars with '_'.
def uniquify_headers(input_fastq_path, out_fastq_path, fp_out):
	try:
		fp_in = open(input_fastq_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
		exit(0);
	num_seqs = 0;
	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		if (len(read) == 0):
			break;
		read[0] = read[0].split()[0] + '-%d' % (num_seqs);
		fp_out.write('\n'.join(read) + '\n');
		num_seqs += 1;
	sys.stderr.write('\n');
	fp_in.close();
	fp_out.close();
	# return header_hash;

if __name__ == "__main__":
	if (len(sys.argv) < 2):
		sys.stderr.write('Various filtering methods for FASTA/FASTQ files.\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\theader\n');
		sys.stderr.write('\tduplicateid\n');
		sys.stderr.write('\tjoin\n');
		sys.stderr.write('\tlen\n');
		sys.stderr.write('\tqualstats\n');
		sys.stderr.write('\treverse\n');
		sys.stderr.write('\thardclip\n');
		sys.stderr.write('\tspecialchars\n');
		sys.stderr.write('\tmarginalign\n');
		sys.stderr.write('\twrap\n');
		sys.stderr.write('\ttouppercase\n');
		sys.stderr.write('\tenumerate\n');
		sys.stderr.write('\tgiheader\n');
		sys.stderr.write('\tuniquify\n');

		exit(0);

	if (sys.argv[1] == 'header'):
		if (len(sys.argv) < 4 or len(sys.argv) > 5):
			sys.stderr.write('Takes a FASTA/FASTQ file and a file containing header patterns. Extracts all sequences containing the specified header patterns.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_pattern_file> <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			sys.stderr.write('\n');
			sys.stderr.write('\tIf any of the patterns match a substring of a header, the sequence will be outputted.\n');
			sys.stderr.write('\tTo match the entire header and not only a substring, place a semicolon ";" character at the end of the pattern line.\n');
			sys.stderr.write('\n');
			exit(0);

		header_patterns_path = sys.argv[2];
		input_fastq_path = sys.argv[3];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 5):
			out_fastq_path = sys.argv[4];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(str(e));
				exit(0);

		filter_seqs_by_header(input_fastq_path, header_patterns_path, out_fastq_path, fp_out);

		if (fp_out != sys.stdout):
			fp_out.close();
		exit(0);

	elif (sys.argv[1] == 'duplicateid'):
		if (len(sys.argv) < 3 or len(sys.argv) > 4):
			sys.stderr.write('Filters the duplicate FASTA/FASTQ entries. Duplicates are equal entries up to the first space. (E.g. gi|49175990|ref|NC_000913.2|)\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			sys.stderr.write('\n');
			exit(0);

		input_fastq_path = sys.argv[2];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 4):
			out_fastq_path = sys.argv[3];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(str(e));
				exit(0);

		filter_duplicate_ncbi_id(input_fastq_path, out_fastq_path, fp_out);

		if (fp_out != sys.stdout):
			fp_out.close();
		exit(0);

	elif (sys.argv[1] == 'join'):
		if (len(sys.argv) < 3 or len(sys.argv) > 4):
			sys.stderr.write('If a FASTA/FASTQ entries are stored in lines of specific width, this tool joins them all into one line.');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			sys.stderr.write('\n');
			exit(0);

		input_fastq_path = sys.argv[2];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 4):
			out_fastq_path = sys.argv[3];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(e);
				exit(0);

		join_fastq_lines(input_fastq_path, out_fastq_path, fp_out);

		if (fp_out != sys.stdout):
			fp_out.close();
		exit(0);

	elif (sys.argv[1] == 'len'):
		if (len(sys.argv) < 4 or len(sys.argv) > 5):
			sys.stderr.write('Extracts only sequences which have length >= min_len.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s min_len <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			sys.stderr.write('\n');
			exit(0);

		min_length = int(sys.argv[2]);
		input_fastq_path = sys.argv[3];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 5):
			out_fastq_path = sys.argv[4];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(e);
				exit(0);

		filter_seqs_by_length(input_fastq_path, min_length, out_fastq_path, fp_out);

		if (fp_out != sys.stdout):
			fp_out.close();
		exit(0);

	elif (sys.argv[1] == 'qualstats'):
		if (len(sys.argv) < 3 or len(sys.argv) > 3):
			sys.stderr.write('This is not an actual filter, but calculates statistics of the base qualities provided in the FASTQ file.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			sys.stderr.write('\n');
			exit(0);

		input_fastq_path = sys.argv[2];

		base_quality_stats(input_fastq_path);

		exit(0);

	elif (sys.argv[1] == 'reverse'):
		if (len(sys.argv) < 3 or len(sys.argv) > 3):
			sys.stderr.write('Reverse complements all sequences in a given file.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			sys.stderr.write('\n');
			exit(0);

		input_fastq_path = sys.argv[2];

		reverse_complement_seqs(input_fastq_path);

		exit(0);

	elif (sys.argv[1] == 'hardclip'):
		if (len(sys.argv) < 4 or len(sys.argv) > 4):
			sys.stderr.write('Trim each sequence from the left and the right end by a given amount of bases.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s num_bases <input_fastq_file>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			sys.stderr.write('\n');
			exit(0);

		num_bases = int(sys.argv[2]);
		input_fastq_path = sys.argv[3];

		hard_clip(input_fastq_path, num_bases);

		exit(0);

	elif (sys.argv[1] == 'specialchars'):
		if (len(sys.argv) < 3 or len(sys.argv) > 4):
			sys.stderr.write('Removes any non-alnum character from the sequences\' header and replaces it with \'_\'.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			exit(0);

		input_fastq_path = sys.argv[2];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 4):
			out_fastq_path = sys.argv[3];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(str(e));
				exit(0);

		remove_special_chars_from_headers(input_fastq_path, out_fastq_path, fp_out);

		if (fp_out != sys.stdout):
			fp_out.close();

		exit(0);

	elif (sys.argv[1] == 'marginalign'):
		if (len(sys.argv) < 3 or len(sys.argv) > 4):
			sys.stderr.write('Removes any non-alnum character from the sequences\' header and replaces it with \'_\'. Also, removes all reads larger than 50kbp.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			exit(0);

		input_fastq_path = sys.argv[2];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 4):
			out_fastq_path = sys.argv[3];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(str(e));
				exit(0);

		filter_for_marginalign(input_fastq_path, out_fastq_path, fp_out);

		if (fp_out != sys.stdout):
			fp_out.close();

		exit(0);

	elif (sys.argv[1] == 'wrap'):
		if (len(sys.argv) < 4 or len(sys.argv) > 5):
			sys.stderr.write('Wraps FASTA/FASTQ files to have line length up to a given number of characters.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s line_length <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			exit(0);

		wrap_length = int(sys.argv[2]);
		input_fastq_path = sys.argv[3];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 5):
			out_fastq_path = sys.argv[4];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(str(e));
				exit(0);

		wrap_fastq_file(input_fastq_path, wrap_length, out_fastq_path, fp_out)

		if (fp_out != sys.stdout):
			fp_out.close();

		exit(0);

	elif (sys.argv[1] == 'touppercase'):
		if (len(sys.argv) < 4 or len(sys.argv) > 4):
			sys.stderr.write('Converts lowercase bases (e.g. masking regions) to upper case.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			exit(0);

		input_fastq_path = sys.argv[2];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 4):
			out_fastq_path = sys.argv[3];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(str(e));
				exit(0);

		convert_to_uppercase(input_fastq_path, out_fastq_path, fp_out)

		if (fp_out != sys.stdout):
			fp_out.close();

		exit(0);

	elif (sys.argv[1] == 'enumerate'):
		if (len(sys.argv) < 4 or len(sys.argv) > 4):
			sys.stderr.write('Replace headers of FASTQ sequences with their ordinal numbers.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			exit(0);

		input_fastq_path = sys.argv[2];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 4):
			out_fastq_path = sys.argv[3];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(str(e));
				exit(0);

		enumerate_headers(input_fastq_path, out_fastq_path, fp_out)

		if (fp_out != sys.stdout):
			fp_out.close();

		exit(0);

	elif (sys.argv[1] == 'giheader'):
		if (len(sys.argv) < 4 or len(sys.argv) > 4):
			sys.stderr.write('If a FASTQ file contains whitespaces, remove everything after the first whitespace.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			exit(0);

		input_fastq_path = sys.argv[2];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 4):
			out_fastq_path = sys.argv[3];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(str(e));
				exit(0);

		keep_gi_header(input_fastq_path, out_fastq_path, fp_out)

		if (fp_out != sys.stdout):
			fp_out.close();

		exit(0);

	elif (sys.argv[1] == 'uniquify'):
		if (len(sys.argv) < 4 or len(sys.argv) > 4):
			sys.stderr.write('Makes all headers in a FASTQ file unique by adding a number at the end of the header. Also, splits headers up to the first whitespace.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
			exit(0);

		input_fastq_path = sys.argv[2];

		out_fastq_path = '';
		fp_out = sys.stdout;
		if (len(sys.argv) == 4):
			out_fastq_path = sys.argv[3];
			if (input_fastq_path == out_fastq_path):
				sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
				exit(0);
			try:
				fp_out = open(out_fastq_path, 'w');
			except Exception, e:
				sys.stderr.write(str(e));
				exit(0);

		uniquify_headers(input_fastq_path, out_fastq_path, fp_out)

		if (fp_out != sys.stdout):
			fp_out.close();

		exit(0);

	else:
		sys.stderr.write('ERROR: Unknown subcommand!\n');
		exit(0);
