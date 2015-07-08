#! /usr/bin/python

# Copyright Ivan Sovic, 2015. www.sovic.org
# 
# Analyzes error rates from a given SAM file, and plots results.

import os;
import sys;
import fastqparser;
import utility_sam;

import numpy as np;

USE_MATPLOTLIB = True;
try:
	import matplotlib.pyplot as plt;
	from matplotlib.font_manager import FontProperties;
except:
	USE_MATPLOTLIB = False;

HIGH_DPI_PLOT = False;
# HIGH_DPI_PLOT = True;

# CIGAR_OPERATIONS = ['M', 'I', 'D', '='];
CIGAR_OPERATIONS_ALL = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'];
CIGAR_OPERATIONS_BASIC = ['M', 'I', 'D', 'S', 'H'];
CIGAR_OPERATIONS_EXTENDED = ['M', 'I', 'D', 'S', 'H', '=', 'X'];

PARAMS_FOR_STATS = ['read_length', 'edit_distance'];
LONG_NAME = {'M': '(Mis)Match', 'I': 'Insertion', 'D': 'Deletion', '=': 'Match', 'X': 'Mismatch', 'read_length': 'Read length', 'edit_distance': 'Edit distance'};



#def SplitCigar(cigar):
	#i = 0;
	#cigarcount_string = '';
	#cigar_operations = [];
	
	#while i < len(cigar):
		#if (cigar[i] in CIGAR_OPERATIONS_EXTENDED):
			#cigar_operations.append([int(cigarcount_string), cigar[i]]);
			#cigarcount_string = '';
		#else:
			#cigarcount_string += cigar[i];
		
		#i += 1;
	
	#return cigar_operations;

def CountCigarOperations(references, sam_line, count_indels_as_events=False):
	if (sam_line.IsMapped() == False):
		return None;
	
	qname = sam_line.qname;
	rname = sam_line.rname;
	cigar = sam_line.cigar;
	read_sequence = sam_line.seq;
	position_reference = sam_line.pos - 1;
	position_read = 0;
	
	ret = None;
	reference = None;
	
	insertions = 0;
	deletions = 0;
	matches = 0;
	mismatches = 0;
	
	try:
		reference = references[rname];			# reference[0] is the header, reference[1] is the sequence.
		#reference_header = reference[0];
		reference_sequence = reference;
	except:
		sys.stderr.write('ERROR: Reference name "%s" not found! Read name: "%s".\n' % (rname, qname));
		return None;
	
	#split_cigar = SplitCigar(cigar);
	split_cigar = sam_line.SplitCigar();
	
	#print reference_header;
	#print ' ';
	
	clipped_read_length = 0;
	
	i = 0;
	while i < len(split_cigar):
		cigarop = split_cigar[i];
		
		#print '[%d] ' % i, cigarop;
		
		if (cigarop[1] == 'H'):
			i += 1;
			continue;
		
		if (cigarop[1] == 'S'):
			position_read += cigarop[0];
			i += 1;
			continue;
		
		if (cigarop[1] == 'I'):
			if (count_indels_as_events == True):
				insertions += 1;
			else:
				insertions += cigarop[0];
			position_read += cigarop[0];
			clipped_read_length += cigarop[0];
			i += 1;
			continue;
		
		if (cigarop[1] == 'D'):
			if (count_indels_as_events == True):
				deletions += 1;
			else:
				deletions += cigarop[0];
			position_reference += cigarop[0];
			i += 1;
			continue;
		
		if (cigarop[1] in ['M', '=', 'X']):
			j = 0;
			while j < cigarop[0]:
				#print 'ref[%d] %s %s %s read[%d]' % (position_reference, reference_sequence[position_reference], cigarop[1], read_sequence[position_read], position_read);
				if reference_sequence[position_reference] == read_sequence[position_read]:
					matches += 1;
				else:
					mismatches += 1;
					
				position_read += 1;
				position_reference += 1;
				clipped_read_length += 1;
				
				j += 1;
			
			
		
		#if (cigarop[1] in ['M', 'I', 'S', '=', 'X']):
			#position += 1;
		
		#if (cigarop[1] in 'IS'):
			#position +=
		
		i += 1;
	
	#print (matches + mismatches + insertions);
	#print (mismatches + insertions + deletions);
	
	read_length = len(read_sequence);
	errors = mismatches + insertions + deletions;
	error_rate = float(errors) / float(clipped_read_length);
	match_rate = float(matches) / float(clipped_read_length);
	mismatch_rate = float(mismatches) / float(clipped_read_length);
	insertion_rate = float(insertions) / float(clipped_read_length);
	deletion_rate = float(deletions) / float(clipped_read_length);

	# if (match_rate < 0.50):
	# 	sys.stderr.write('\n' + sam_line.FormatAccuracy() + '\n');
	# 	sys.stderr.write('matches = %d, clipped_read_length = %d\n' % (matches, clipped_read_length));

	# print 'mismatches = ', mismatches, ', insertions = ', insertions, ', deletions = ', deletions, ', matches = ', matches, ', errors = ', errors;
	# print '[match_rate, mismatch_rate, insertion_rate, deletion_rate, error_rate, matches, mismatches, insertions, deletions, errors, read_length, clipped_read_length] = ', [match_rate, mismatch_rate, insertion_rate, deletion_rate, error_rate, matches, mismatches, insertions, deletions, errors, read_length, clipped_read_length];

#    if (error_rate > 0.43):
#            sys.stdout.write('%s\t%d\n' % (sam_line.qname, sam_line.clipped_pos));
#            sys.stdout.flush();

	# if (error_rate > 0.49 and error_rate < 0.51):
	# 	print 'errors = ', errors;
	# 	sam_line.Verbose();
	# 	sys.stdout.flush();
	# 	exit(1);
	
	#exit(1);
	
	return [match_rate, mismatch_rate, insertion_rate, deletion_rate, error_rate, matches, mismatches, insertions, deletions, errors, read_length, clipped_read_length];

def FilterZeros(hist_values):
	ret_hist_x = [];
	ret_hist_values = [];
	
	if (len(hist_values) > 0):
		ret_hist_x.append(0);
		ret_hist_values.append(hist_values[0]);
		
	i = 1;
	while i < (len(hist_values) - 1):
		if (hist_values[i] == 0 and hist_values[i-1] != 0 and hist_values[i+1] != 0):
			i += 1;
			continue;
		ret_hist_x.append(i);
		ret_hist_values.append(hist_values[i]);
		i += 1;
	
	if (len(hist_values) > 1):
		ret_hist_x.append((len(hist_values) - 1));
		ret_hist_values.append(hist_values[len(hist_values) - 1]);
	
	#ret_hist_x = range(0, len(hist_values));
	#ret_hist_values = hist_values;
	
	return [ret_hist_x, ret_hist_values];
	
	

def PlotErrorRates(error_rate_hist, insertion_hist, deletion_hist, snp_hist, match_hist, dataset_name, out_png_path):
	fig = plt.figure();
	ax1 = plt.subplot(111);
	
	[xvalues_error_rate, y_error_rate_hist] = FilterZeros(error_rate_hist);
	[xvalues_insertion, y_insertion_hist] = FilterZeros(insertion_hist);
	[xvalues_deletion, y_deletion_hist] = FilterZeros(deletion_hist);
	[xvalues_snp, y_snp_hist] = FilterZeros(snp_hist);
	[xvalues_match_rate, y_match_hist] = FilterZeros(match_hist);
	
	
	
	xvalues = range(0, 101);
	ax1.plot(xvalues_insertion, y_insertion_hist, label='Insertions (I)', linewidth=2);
	ax1.plot(xvalues_deletion, y_deletion_hist, label='Deletions (D)', linewidth=2);
	ax1.plot(xvalues_snp, y_snp_hist, label='Mismatches (M)', linewidth=2);
	ax1.plot(xvalues_error_rate, y_error_rate_hist, label='I + D + M', linewidth=2);
	ax1.plot(xvalues_match_rate, y_match_hist, 'k--', label='Matches', linewidth=2);
	
	fontP = FontProperties()
	# fontP.set_size('small')
	fontP.set_size(14)
	ax1.grid();
	ax1.legend(prop=fontP, loc='upper right');

#	title_string = 'Analysis of error rates as obtained with GraphMap\n%s' % dataset_name;
	# modified_dataset_name = dataset_name;
	# if ('graphmap' in dataset_name.lower()):
	# 	modified_dataset_name = 'GraphMap';
	# elif ('last' in dataset_name.lower()):
	# 	modified_dataset_name = 'LAST';
	# elif ('bwamem' in dataset_name.lower()):
	# 	modified_dataset_name = 'BWA-MEM';
	# elif ('blasr' in dataset_name.lower()):
	# 	modified_dataset_name = 'BLASR';
	# title_string = 'Analysis of error rates\nReads mapped with %s' % modified_dataset_name;
	# ax1.set_title(title_string);
	
	ax1.set_xlabel('Error rate [%]');
	ax1.set_ylabel('Alignments [%]');

	font = {'family' : 'sans-serif',
		'weight' : 'normal',
		'size'   : 18}
	plt.rc('font', **font)


	# plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);
	if (HIGH_DPI_PLOT == False):
		plt.savefig(out_png_path, bbox_inches='tight'); # , dpi=1000);
	else:
		plt.savefig(out_png_path, bbox_inches='tight', dpi=1000);


	
	ax1.grid();
	ax1.legend(prop=fontP, loc='lower right');
	ax1.grid();

def CollectAccuracy(sam_path, accuracy_path, suppress_error_messages=False):
	sam_basename = os.path.splitext(os.path.basename(sam_path))[0];
	out_png_path = accuracy_path + '.png';

	error_rate_hist = [0 for i in range(0, 101)];
	insertion_hist = [0 for i in range(0, 101)];
	deletion_hist = [0 for i in range(0, 101)];
	snp_hist = [0 for i in range(0, 101)];
	match_hist = [0 for i in range(0, 101)];
	
	average_error_rate = 0.0;		std_error_rate = 0.0;
	average_insertion_rate = 0.0;		std_insertion_Rate = 0.0;
	average_deletion_rate = 0.0;		std_deletion_rate = 0.0;
	average_snp_rate = 0.0;			std_snp_rate = 0.0;
	average_match_rate = 0.0;		std_match_rate = 0.0;
	
	all_error_rates = [];
	all_insertion_rates = [];
	all_deletion_rates = [];
	all_mismatch_rates = [];
	all_match_rates = [];
	all_read_lengths = [];
	
	fp = None;
	try:
		fp = open(accuracy_path, 'r');
	except IOError:
		if (suppress_error_messages == False):
			sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!\n' % (__name__, accuracy_path));
		return ['', error_rate_hist, insertion_hist, deletion_hist, snp_hist, match_hist, sam_basename, out_png_path];
	
	num_lines = 0;
	for line in fp:
		line = line.rstrip();
		if (len(line) == 0):
			continue;
		
		split_line = line.split('\t');
		#print split_line;
		
		#[match_rate, mismatch_rate, insertion_rate, deletion_rate, error_rate, matches, mismatches, insertions, deletions, errors, read_length, clipped_read_length] = split_line;
		error_rate = float(split_line[4]);
		insertion_rate = float(split_line[2]);
		deletion_rate = float(split_line[3]);
		mismatch_rate = float(split_line[1]);
		match_rate = float(split_line[0]);
		read_length = float(split_line[10]);
		
		if (error_rate >= 1.0 or insertion_rate >= 1.0 or deletion_rate >= 1.0 or mismatch_rate >= 1.0 or match_rate >= 1.0):
			continue;
		
		#print error_rate
		error_rate_hist[int(error_rate * 100)] += 1;
		insertion_hist[int(insertion_rate * 100)] += 1;
		deletion_hist[int(deletion_rate * 100)] += 1;
		snp_hist[int(mismatch_rate * 100)] += 1;
		match_hist[int(match_rate * 100)] += 1;
		
		all_error_rates.append(error_rate);
		all_insertion_rates.append(insertion_rate);
		all_deletion_rates.append(deletion_rate);
		all_mismatch_rates.append(mismatch_rate);
		all_match_rates.append(match_rate);
		all_read_lengths.append(read_length);
		
		#average_error_rate += error_rate;
		#average_insertion_rate += insertion_rate;
		#average_deletion_rate += deletion_rate;
		#average_snp_rate += mismatch_rate;
		#average_match_rate += match_rate;
		
		#if (error_rate == 0.0):
			#print 'Alignment %d has error_rate equal to 0.0!' % num_lines;
		num_lines += 1;
	
	fp.close();
	
	column_labels = ['mean', 'std', 'median', 'min', 'max'];
	
	error_rate_stats = [np.mean(all_error_rates), np.std(all_error_rates), np.median(all_error_rates), np.min(all_error_rates), np.max(all_error_rates)];
	insertion_rate_stats = [np.mean(all_insertion_rates), np.std(all_insertion_rates), np.median(all_insertion_rates), np.min(all_insertion_rates), np.max(all_insertion_rates)];
	deletion_rate_stats = [np.mean(all_deletion_rates), np.std(all_deletion_rates), np.median(all_deletion_rates), np.min(all_deletion_rates), np.max(all_deletion_rates)];
	mismatch_rate_stats = [np.mean(all_mismatch_rates), np.std(all_mismatch_rates), np.median(all_mismatch_rates), np.min(all_mismatch_rates), np.max(all_mismatch_rates)];
	match_rate_stats = [np.mean(all_match_rates), np.std(all_match_rates), np.median(all_match_rates), np.min(all_match_rates), np.max(all_match_rates)];
	read_length_stats = [np.mean(all_read_lengths), np.std(all_read_lengths), np.median(all_read_lengths), np.min(all_read_lengths), np.max(all_read_lengths)];
	
	# ret_lines = '';
	# ret_lines += 'Error rate stats:     \t%s\n' % (',\t'.join(['%s = %.2f' % (column_labels[i], error_rate_stats[i]) for i in range(len(error_rate_stats))]));
	# ret_lines += 'Insertion rate stats: \t%s\n' % (',\t'.join(['%s = %.2f' % (column_labels[i], insertion_rate_stats[i]) for i in range(len(error_rate_stats))]));
	# ret_lines += 'Deletion rate stats:  \t%s\n' % (',\t'.join(['%s = %.2f' % (column_labels[i], deletion_rate_stats[i]) for i in range(len(error_rate_stats))]));
	# ret_lines += 'Mismatch rate stats:  \t%s\n' % (',\t'.join(['%s = %.2f' % (column_labels[i], mismatch_rate_stats[i]) for i in range(len(error_rate_stats))]));
	# ret_lines += 'Match rate stats:     \t%s\n' % (',\t'.join(['%s = %.2f' % (column_labels[i], match_rate_stats[i]) for i in range(len(error_rate_stats))]));
	# ret_lines += 'Read length stats:    \t%s\n' % (',\t'.join(['%s = %.2f' % (column_labels[i], read_length_stats[i]) for i in range(len(error_rate_stats))]));

	
	ret_lines = '';
	lines_error_rates = '';
	lines_error_rates += '                      \t%s\n' % ('\t'.join(['%s' % (column_labels[i]) for i in range(len(error_rate_stats))]));
	lines_error_rates += 'Error rate:     \t%s\n' % ('\t'.join(['%.2f' % (error_rate_stats[i]) for i in range(len(error_rate_stats))]));
	lines_error_rates += 'Insertion rate: \t%s\n' % ('\t'.join(['%.2f' % (insertion_rate_stats[i]) for i in range(len(error_rate_stats))]));
	lines_error_rates += 'Deletion rate:  \t%s\n' % ('\t'.join(['%.2f' % (deletion_rate_stats[i]) for i in range(len(error_rate_stats))]));
	lines_error_rates += 'Mismatch rate:  \t%s\n' % ('\t'.join(['%.2f' % (mismatch_rate_stats[i]) for i in range(len(error_rate_stats))]));
	lines_error_rates += 'Match rate:     \t%s\n' % ('\t'.join(['%.2f' % (match_rate_stats[i]) for i in range(len(error_rate_stats))]));
	lines_error_rates += 'Read length:    \t%s\n' % ('\t'.join(['%.2f' % (read_length_stats[i]) for i in range(len(error_rate_stats))]));

	ret_lines += lines_error_rates;
	
	#ret_lines += 'Error rate stats: %s\n' % (', '.join(['%.2f' % value for value in error_rate_stats]));
	#ret_lines += 'Insertion rate stats: %s\n' % (', '.join(['%.2f' % value for value in insertion_rate_stats]));
	#ret_lines += 'Deletion rate stats: %s\n' % (', '.join(['%.2f' % value for value in deletion_rate_stats]));
	#ret_lines += 'Mismatch rate stats: %s\n' % (', '.join(['%.2f' % value for value in mismatch_rate_stats]));
	#ret_lines += 'Match rate stats: %s\n' % (', '.join(['%.2f' % value for value in match_rate_stats]));
	#ret_lines += 'Read length stats: %s\n' % (', '.join(['%.2f' % value for value in read_length_stats]));

	insertion_ratio = int((insertion_rate_stats[0]/error_rate_stats[0]) * 100);
	deletion_ratio = int((deletion_rate_stats[0]/error_rate_stats[0]) * 100);
	mismatch_ratio = 100 - insertion_ratio - deletion_ratio;
	
	ret_lines += 'Difference ratio: %d:%d:%d (mismatch:insertion:deletion)\n' % (mismatch_ratio, insertion_ratio, deletion_ratio);
	
	#print ret_lines;

	#if (num_lines > 0):
		#average_error_rate /= num_lines;
		#average_insertion_rate /= num_lines;
		#average_deletion_rate /= num_lines;
		#average_snp_rate /= num_lines;
		#average_match_rate /= num_lines;
	
	#i = 0;
	#while i < (len(error_rate_hist)):
		#if (error_rate_hist[i] < 10):
			#print 'error_rate_hist[%d] = %d' % (i, error_rate_hist[i]);
		#i += 1;
	
	out_hist_path = accuracy_path + '.plot';
	try:
		fp_hist = open(out_hist_path, 'w');
		fp_hist.write(sam_path + '\n');
		fp_hist.write(accuracy_path + '\n');
		fp_hist.write('x\t-\terror_rate\t' + '\t'.join([str(value) for value in range(0, 101)]) + '\n');
		fp_hist.write('y\t' + sam_basename + '\terror_rate_count\t' + '\t'.join([str(value) for value in error_rate_hist]) + '\n');
		fp_hist.write('y\t' + sam_basename + '\tinsertion_rate_count\t' + '\t'.join([str(value) for value in insertion_hist]) + '\n');
		fp_hist.write('y\t' + sam_basename + '\tdeletion_rate_count\t' + '\t'.join([str(value) for value in deletion_hist]) + '\n');
		fp_hist.write('y\t' + sam_basename + '\tsnp_rate_count\t' + '\t'.join([str(value) for value in snp_hist]) + '\n');
		fp_hist.close();
	except IOError:
		if (suppress_error_messages == False):
			sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % out_hist_path);

	### Normalizing the histograms.
	sum_error_rate_counts = sum(error_rate_hist);
	error_rate_hist = [(100*float(value)/float(sum_error_rate_counts)) for value in error_rate_hist];

	sum_insertion_rate_counts = sum(insertion_hist);
	insertion_hist = [(100*float(value)/float(sum_insertion_rate_counts)) for value in insertion_hist];

	sum_deletion_rate_counts = sum(deletion_hist);
	deletion_hist = [(100*float(value)/float(sum_deletion_rate_counts)) for value in deletion_hist];

	sum_snp_rate_counts = sum(snp_hist);
	snp_hist = [(100*float(value)/float(sum_snp_rate_counts)) for value in snp_hist];

	sum_match_rate_counts = sum(match_hist);
	match_hist = [(100*float(value)/float(sum_match_rate_counts)) for value in match_hist];



	PlotErrorRates(error_rate_hist, insertion_hist, deletion_hist, snp_hist, match_hist, sam_basename, out_png_path);
	
	#return [match_rate, mismatch_rate, insertion_rate, deletion_rate, error_rate, matches, mismatches, insertions, deletions, errors, read_length, clipped_read_length];
	
	return [ret_lines, lines_error_rates, error_rate_hist, insertion_hist, deletion_hist, snp_hist, match_hist, sam_basename, out_png_path];
	
	

def ProcessSAM(references, sam_path, accuracy_counts_path, count_indels_as_events=False):
	i = 0;
			
	accuracy_counts = [];

	fp_sam = None;
	fp_accuracy_counts = None;

	try:
                if sam_path == "-":
                        fp_sam = sys.stdin
                else:
                        fp_sam = open(sam_path, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % sam_path);
		exit(1);

	try:
		fp_accuracy_counts = open(accuracy_counts_path, 'w');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % accuracy_counts_path);
                if fp_sam != sys.stdin:
                        fp_sam.close();
		exit(1);
		
	for line in fp_sam:
		i += 1;
		
		if ((i % 100) == 0):
			sys.stderr.write('Processing alignment: %d...\r' % (i));
			sys.stderr.flush();

		line = line.rstrip();
		
		if (len(line) == 0 or line[0] == '@'):
			continue;
		
		sam_line = utility_sam.SAMLine(line);
		
		single_counts = CountCigarOperations(references, sam_line, count_indels_as_events);
		
		if (single_counts != None):
			accuracy_counts.append(single_counts);
			fp_accuracy_counts.write('\t'.join([('%.2f' % float(value)) for value in single_counts]) + '\t' + sam_line.qname + '\n');

	sys.stderr.write('\n');
	sys.stderr.flush();

        if fp_sam != sys.stdin:
                fp_sam.close();
	fp_accuracy_counts.close();

def ProcessFromFiles(reference_file, sam_path, out_accuracy_counts_path, count_indels_as_events=False):
	[ref_headers, ref_seqs, ref_quals] = fastqparser.read_fastq(reference_file);
	references = {};
	accuracy_counts = [];
	
	i = 0;
	while i < len(ref_headers):
		header = ref_headers[i];
		seq = ref_seqs[i];
		references[header] = seq;
		references[header.split()[0]] = seq;
		i += 1;
	ProcessSAM(references, sam_path, out_accuracy_counts_path, count_indels_as_events);

def print_usage_and_exit():
	sys.stderr.write('Analyze error rates from a given SAM file.\n');
	sys.stderr.write('Usage:\n');
	sys.stderr.write('\t%s mode <reference_fasta> <input_sam_file>\n' % (sys.argv[0]));
	sys.stderr.write('\n');
	sys.stderr.write('mode - "base" calculates error rates with indels viewed as per-base events.\n');
	sys.stderr.write('       "event" calculates error rates with indels counted as entire events.\n');
	sys.stderr.write('       "collect" does not recalculate error rates, but just collects and plots.\n');
	sys.stderr.write('\n');
	exit(0);

if __name__ == "__main__":
	if (len(sys.argv) != 4):
		print_usage_and_exit();

	mode = sys.argv[1];
	reference_file = sys.argv[2];
	sam_file = sys.argv[3];

	if (not mode in ['base', 'event', 'collect']):
		print_usage_and_exit();

	out_path = os.path.dirname(sam_file);
	if (out_path == ''):
		out_path = './';

	# sys.stderr.write('Using mode: %s\n' % mode);
	
	out_accuracy_counts_path = '%s/%s-error_rates-%s.csv' % (out_path, os.path.basename(sam_file), mode);
	# out_summary_path = '%s/%s-error_rates-%s-summary.txt' % (out_path, os.path.basename(sam_file), mode);
	out_summary_stats_path = '%s/%s-error_rates-%s-summary.csv' % (out_path, os.path.basename(sam_file), mode);

	if (mode == 'base'):
		ProcessFromFiles(reference_file, sam_file, out_accuracy_counts_path, False);
	if (mode == 'event'):
		ProcessFromFiles(reference_file, sam_file, out_accuracy_counts_path, True);

	[ret_lines, lines_error_rates, error_rate_hist, insertion_hist, deletion_hist, snp_hist, match_hist, sam_basename, out_png_path] = CollectAccuracy(sam_file, out_accuracy_counts_path);
	sys.stdout.write(ret_lines);

	# fp_summary = open(out_summary_path, 'w');
	# fp_summary.write(ret_lines);
	# fp_summary.close();

	fp_summary_stats_path = open(out_summary_stats_path, 'w');
	fp_summary_stats_path.write(lines_error_rates);
	fp_summary_stats_path.close();

	plt.show();
