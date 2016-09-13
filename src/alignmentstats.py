#! /usr/bin/env python

import os;
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
#GOLDEN_PATH = SCRIPT_PATH + '/../../../golden-bundle';
import sys;
# sys.path.append(SCRIPT_PATH + '/scripts')
#sys.path.append(GOLDEN_PATH + '/src')
import subprocess;
import numpy as np;

# import filesandfolders;
import fastqparser;
import utility_sam;
import errorrates;
import count_mapped_reads;
# import consensus_stats;
import consensus;
# import memtime_stats;

USE_MATPLOTLIB = True;
try:
        import matplotlib;
        matplotlib.use('Agg');
        import matplotlib.pyplot as plt;
        from matplotlib.font_manager import FontProperties;
except Exception, e:
        USE_MATPLOTLIB = False;
        sys.stderr.write('Exception when importing Matplotlib.\n');
        sys.stderr.write(str(e) + '\n');
        sys.stderr.write('\n');

COLLECT_RESULTS_WITHOUT_REANALYSIS = False;
#COLLECT_RESULTS_WITHOUT_REANALYSIS = True;

MODE_CODE_CALC_CONSENSUS = (1 << 0);
MODE_CODE_CALC_ERROR_RATE = (1 << 1);
MODE_CODE_CALC_ERROR_RATE_INDELS_AS_EVENTS = (1 << 2);
MODE_CODE_CALC_READ_COUNTS = (1 << 3);
MODE_CODE_COLLECT_CONSENSUS = (1 << 4);
MODE_CODE_COLLECT_ERROR_RATE = (1 << 5);
MODE_CODE_COLLECT_ERROR_RATE_INDELS_AS_EVENTS = (1 << 6);
MODE_CODE_COLLECT_READ_COUNTS = (1 << 7);
MODE_CODE_HEADER = (1 << 8);

def print_usage_and_exit():
	print 'Usage:';
	print '\t%s file|folder mode_code <sam_path> <reference_path> <reads_path> [consensus_coverage_threshold] [<simulated_sam_path>]' % sys.argv[0];
	print '';
	print '\tmode_code\t0 Run all analyses.';
	print '\t\t\t1 Only collect the results';
	print '\t\t\thcalc';
	print '\t\t\tcalc';
	print '\t\t\thcollect';
	print '\t\t\tcollect';
	print '\t\t\thcalccons';
	print '\t\t\tcalccons';
	print '\t\t\thcalcerrors';
	print '\t\t\tcalcerrors';
	print '\t\t\t';
	exit(1);

def run_from_args(cmd_args):
	file_or_folder = cmd_args[0];
	mode = cmd_args[1];
	sam_path = cmd_args[2];
	reference_file = cmd_args[3];
	reads_path = cmd_args[4];
	simulated_sam_path = '';

	consensus_coverage_threshold = 20;
	if (len(sys.argv) >=7 ):
		consensus_coverage_threshold = int(sys.argv[6]);

	if (len(cmd_args) >= 6):
		simulated_sam_path = int(cmd_args[5]);
	if (len(cmd_args) >= 7):
		simulated_sam_path = cmd_args[6];

	if (mode == 'calc'):
		mode_code = MODE_CODE_CALC_CONSENSUS | MODE_CODE_CALC_ERROR_RATE | MODE_CODE_CALC_ERROR_RATE_INDELS_AS_EVENTS | MODE_CODE_CALC_READ_COUNTS;
	elif (mode == 'hcalc'):
		mode_code = MODE_CODE_CALC_CONSENSUS | MODE_CODE_CALC_ERROR_RATE | MODE_CODE_CALC_ERROR_RATE_INDELS_AS_EVENTS | MODE_CODE_CALC_READ_COUNTS | MODE_CODE_HEADER;
	elif (mode == 'collect'):
		mode_code = MODE_CODE_COLLECT_CONSENSUS | MODE_CODE_COLLECT_ERROR_RATE | MODE_CODE_COLLECT_ERROR_RATE_INDELS_AS_EVENTS | MODE_CODE_COLLECT_READ_COUNTS;
	elif (mode == 'hcollect'):
		mode_code = MODE_CODE_COLLECT_CONSENSUS | MODE_CODE_COLLECT_ERROR_RATE | MODE_CODE_COLLECT_ERROR_RATE_INDELS_AS_EVENTS | MODE_CODE_COLLECT_READ_COUNTS | MODE_CODE_HEADER;

	elif (mode == 'hcalccons'):
		mode_code = MODE_CODE_CALC_CONSENSUS | MODE_CODE_COLLECT_ERROR_RATE | MODE_CODE_COLLECT_ERROR_RATE_INDELS_AS_EVENTS | MODE_CODE_COLLECT_READ_COUNTS | MODE_CODE_HEADER;
	elif (mode == 'calccons'):
		mode_code = MODE_CODE_CALC_CONSENSUS | MODE_CODE_COLLECT_ERROR_RATE | MODE_CODE_COLLECT_ERROR_RATE_INDELS_AS_EVENTS | MODE_CODE_COLLECT_READ_COUNTS;

	elif (mode == 'hcalcerrors'):
		mode_code = MODE_CODE_COLLECT_CONSENSUS | MODE_CODE_CALC_ERROR_RATE | MODE_CODE_CALC_ERROR_RATE_INDELS_AS_EVENTS | MODE_CODE_COLLECT_READ_COUNTS | MODE_CODE_HEADER;
	elif (mode == 'calcerrors'):
		mode_code = MODE_CODE_COLLECT_CONSENSUS | MODE_CODE_CALC_ERROR_RATE | MODE_CODE_CALC_ERROR_RATE_INDELS_AS_EVENTS | MODE_CODE_COLLECT_READ_COUNTS;

	else:
		# mode_code = MODE_CODE_CALC_CONSENSUS | MODE_CODE_CALC_ERROR_RATE | MODE_CODE_CALC_READ_COUNTS;
		mode_code = int(mode);

	if (file_or_folder == 'file'):
		# output_folder_path = os.path.dirname(sam_path);
#		analyze_single_sam(sam_path, reference_file, output_folder_path, num_reads_in_reference_sam=1, suppress_error_messages=False, consensus_coverage_thresholds=[20]);
		# analyze_single_sam(mode_code, sam_path, reference_file, reads_path, simulated_sam_path);
		# sys.stdout.write(analyze_single_sam((MODE_CODE_CALC_READ_COUNTS | MODE_CODE_COLLECT_CONSENSUS), sam_file, reference_file, reads_file, simulated_sam_path=simulated_sam_path, consensus_coverage_threshold=consensus_coverage_threshold));
		sys.stdout.write(analyze_single_sam(mode_code, sam_path, reference_file, reads_path, simulated_sam_path=simulated_sam_path, consensus_coverage_threshold=consensus_coverage_threshold));
	elif (file_or_folder == 'folder'):
		# analyze_all_sams_from_folder(sam_path, reference_file);
		sys.stderr.write('Warning: Analysis of all SAMs in a given folder not implemented yet! Skipping.\n');
		pass;
	else:
		print_usage_and_exit();
	
def ParseMemTime(sam_file):
	memtime_path = os.path.splitext(sam_file)[0] + '.memtime';
	cmdline = '';
	realtime = 0.0;
	cputime = 0.0;
	usertime = 0.0;
	systemtime = 0.0;
	maxrss = 0.0;
	rsscache = 0.0;
	lines = [];
	try:
		fp = open(memtime_path, 'r');
		lines = [line.strip() for line in fp.readlines() if (len(line.strip()) > 0)];
		fp.close();

		if (len(lines) > 0):
			cmdline = lines[0];
			realtime = float(lines[1].split(':')[1].strip().split(' ')[0].strip());
			cputime = float(lines[2].split(':')[1].strip().split(' ')[0].strip());
			usertime = float(lines[3].split(':')[1].strip().split(' ')[0].strip());
			systemtime = float(lines[4].split(':')[1].strip().split(' ')[0].strip());
			maxrss = float(lines[5].split(':')[1].strip().split(' ')[0].strip());
			if (len(lines) > 6):
				rsscache = float(lines[6].split(':')[1].strip().split(' ')[0].strip());
			else:
				rsscache = -1.0;
	except IOError:
		sys.stderr.write('ERROR: Could not parse memtime stats from "%s"!\n' % memtime_path);

	return [cmdline, realtime, cputime, usertime, systemtime, maxrss, rsscache, lines];

def ParseVariantStats(variant_summary_path):
	alignments_file = '';
	mpileup_file = '';
	coverage_threshold = 0;
	snp_count = 0;	insertion_count = 0;	deletion_count = 0;
	num_undercovered_bases = 0;	num_called_bases = 0;	num_correct_bases = 0;
	average_coverage = 0.0;
	lines = [];

	try:
		fp = open(variant_summary_path, 'r');
		lines = [line.strip() for line in fp.readlines() if (len(line.strip()) > 0)];
		fp.close();
		for line in lines:
			split_line = line.split(':');
			param_name = split_line[0].strip();
			param_value = split_line[1].strip();
			if (param_name == 'alignments_file'):
				alignments_file = param_value;
			if (param_name == 'mpileup_file'):
				mpileup_file = param_value;
			if (param_name == 'coverage_threshold'):
				coverage_threshold = int(param_value);
			if (param_name == 'snp_count'):
				snp_count = int(param_value);
			if (param_name == 'insertion_count'):
				insertion_count = int(param_value);
			if (param_name == 'deletion_count'):
				deletion_count = int(param_value);
			if (param_name == 'num_undercovered_bases'):
				num_undercovered_bases = int(param_value);
			if (param_name == 'num_called_bases'):
				num_called_bases = int(param_value);
			if (param_name == 'num_correct_bases'):
				num_correct_bases = int(param_value);
			if (param_name == 'average_coverage'):
				average_coverage = float(param_value);
	except IOError:
		sys.stderr.write('ERROR: Could not parse variant stats from "%s"!\n' % variant_summary_path);

	return [alignments_file, mpileup_file, coverage_threshold, snp_count, insertion_count, deletion_count, num_undercovered_bases, num_called_bases, num_correct_bases, average_coverage, lines];

# 12
# 76
# analyze_single_sam((MODE_CODE_CALC_READ_COUNTS | MODE_CODE_COLLECT_CONSENSUS), sam_file, reference_file, reads_file, simulated_sam_path='', consensus_coverage_threshold=20); 

# def Main(sam_file='', reference_file='', output_path='', num_reads_in_reference_sam=1, suppress_error_messages=False, consensus_coverage_thresholds=[20]):
def analyze_single_sam(mode_code, sam_file, reference_file, reads_file, simulated_sam_path='', consensus_coverage_threshold=20):
	
	if (sam_file == '' or reference_file == ''):
		print '[analyze_single_sam] ERROR: No input files given.';
		print_usage_and_exit();
		return;

	dir_name = os.path.dirname(sam_file);
	if (dir_name == ''):
		dir_name = '.';

	# Create the output path if it doesn't exist yet.
	output_folder_intermediate = dir_name + '/analysis-intermediate';
	output_folder_final = dir_name + '/analysis-final';
	if not os.path.exists(output_folder_intermediate):
		os.makedirs(output_folder_intermediate);
	if not os.path.exists(output_folder_final):
		os.makedirs(output_folder_final);
	
	dataset_name = os.path.splitext(os.path.basename(sam_file))[0];
	summary = '';
	csv_line = '%s\t' % (dataset_name);
	csv_header = 'mapper\t';

	# Defining output filenames.
	out_accuracy_counts_path = '%s/error_rates-individualbases-%s.csv' % (output_folder_intermediate, dataset_name);
	out_accuracy_counts_indel_events_path = '%s/error_rates-eventsindel-%s.csv' % (output_folder_intermediate, dataset_name);
	consensus_prefix = '%s/consensus-%s' % (output_folder_intermediate, dataset_name);
	out_summary_path = '%s/summary_sam_analysis-%s.txt' % (output_folder_intermediate, dataset_name);
	out_consensus_plot_lines = consensus_prefix + '.plot';
	out_consensus_plot_png = consensus_prefix + '.plot.png';
	out_count_mapped_reads_prefix = '%s/count_reads-%s' % (output_folder_intermediate, dataset_name);
	out_results_path = '%s/results-%s.csv' % (output_folder_final, dataset_name);


	
	try:
		fp = open(out_summary_path, 'w');
	except IOError:
		sys.stderr.write('ERROR: Could not open summary path "%s" for writing!\n' % out_summary_path);
		os.exit(1);

	# Verbose of the filenames for the output.
	summary_file_paths = '';
	summary_file_paths += 'Analyzing SAM file %s...\n' % (sam_file);
	summary_file_paths += 'Reference file: %s\n' % (reference_file);
	# summary_file_paths += 'Output folder: %s\n' % (output_folder);
	summary_file_paths += 'Accuracy counts: %s\n' % (out_accuracy_counts_path);
	summary_file_paths += 'Consensus prefix: %s\n' % (consensus_prefix);
	summary_file_paths += 'Count mapped reads: %s\n' % (out_count_mapped_reads_prefix);
	summary += '[Paths]\n' + summary_file_paths + '\n';
	# Just simply verbose to screen and summary file.
	sys.stderr.write('[Paths]\n');
	sys.stderr.write(summary_file_paths + '\n');
	fp.write('[Paths]\n');
	fp.write(summary_file_paths + '\n');
	
	# Get the headers from the SAM file (if they exist). Useful if the commandline was stored.
	sam_headers = utility_sam.LoadOnlySAMHeaders(sam_file, False);
	summary_file_sam_headers = '%s\n' % ('\n'.join(sam_headers));
	summary += '[SAM headers]\n' + summary_file_sam_headers + '\n';
	# Just simply verbose to screen and summary file.
	sys.stderr.write('[SAM headers]\n');
	sys.stderr.write(summary_file_sam_headers + '\n');
	fp.write('[SAM headers]\n');
	fp.write(summary_file_sam_headers + '\n');

	# Calculate the input reads statistics.
	[fastqinfo_string, fastqinfo_num_seqs, fastqinfo_total_seq_len, fastqinfo_average_seq_len, temp_max_seq_len] = fastqparser.count_seq_length(reads_file);
	summary_reads_file = '';
	summary_reads_file += 'Number of reads in the input file: %d\n' % (fastqinfo_num_seqs);
	summary_reads_file += 'Total number of bases in the input reads file: %d\n' % (fastqinfo_total_seq_len);
	summary_reads_file += 'Average read length in the input file: %d\n' % (fastqinfo_average_seq_len);
	summary += '[Input reads file]\n' + summary_reads_file + '\n';
	###############################3
	# This version is more verbose, but is different than figure1.xlsx!
	# csv_line += '%d\t%d\t' % (fastqinfo_num_seqs, fastqinfo_total_seq_len);
	# csv_header += 'Num input reads\tNum input bases\t';
	###############################3
	# Just simply verbose to screen and summary file.
	sys.stderr.write('[Input reads file]\n');
	sys.stderr.write(summary_reads_file + '\n');
	fp.write('[Input reads file]\n');
	fp.write(summary_reads_file + '\n');
	
	# header = 'CPU time\tMax RSS\tNum input reads\tNum input bases\t';



	if (mode_code & MODE_CODE_CALC_CONSENSUS):
		summary_consensus = '';
		sys.stderr.write('[main_sam_analysis.py] Calculating consensus statistics, base coverage threshold %d...\n' % (consensus_coverage_threshold));
		consensus.main(sam_file, reference_file, consensus_coverage_threshold, consensus_prefix, thread_id=0);
		sys.stderr.write('\n');
	if ((mode_code & MODE_CODE_CALC_CONSENSUS) or (mode_code & MODE_CODE_COLLECT_CONSENSUS)):
		# Collecting stats from the consensus caller.
		# consensus_plot_lines and consensus_statistics have the same information, only the
		# consensus_statistics object is a dict and the consensus_plot_lines is a formatted string
		# meant for outputting to a summary text file.
		variant_summary_path = ('%s-cov_%d.variant.sum' % (consensus_prefix, consensus_coverage_threshold));
		[alignments_file, mpileup_file, coverage_threshold, snp_count, insertion_count, deletion_count, \
						num_undercovered_bases, num_called_bases, num_correct_bases, average_coverage, variant_lines] = ParseVariantStats(variant_summary_path);

		summary_consensus = '\n'.join(variant_lines);
		summary += '[Consensus statistics]\n' + summary_consensus + '\n\n';

		###############################3
		# This version is more verbose, but is different than figure1.xlsx!
		# csv_line += '%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t' % (	snp_count, deletion_count, insertion_count,
		# 													num_undercovered_bases, num_called_bases,
		# 													num_correct_bases, average_coverage, coverage_threshold);
		# csv_header += 'snp_count\tdeletion_count\tinsertion_count\tnum_undercovered_bases\tnum_called_bases\tnum_correct_bases\taverage_coverage\tcoverage_threshold\t';
		###############################3

		###############################3
		# This version is in the same order as figure1.xlsx!
		csv_line += '%d\t%d\t%d\t%d\t%d\t%d\t%f\t' % (	insertion_count, deletion_count, snp_count,
															num_undercovered_bases, num_called_bases,
															num_correct_bases, average_coverage);
		csv_header += 'Insertions\tDeletions\tSNPs\tUncalled Bases\tnum_called_bases\tnum_correct_bases\taverage_coverage\t';
		###############################3

		# [summary_consensus, consensus_plot_lines, consensus_statistics] = consensus_stats.CollectConsensus(sam_file, consensus_prefix, consensus_coverage_thresholds, suppress_error_messages);
		sys.stderr.write('[Consensus statistics]\n');
		sys.stderr.write(summary_consensus + '\n\n');
		fp.write('[Consensus statistics]\n');
		fp.write(summary_consensus + '\n\n');
		# fp_consensus_plot_lines = open(out_consensus_plot_lines, 'w');
		# fp_consensus_plot_lines.write(consensus_plot_lines);
		# fp_consensus_plot_lines.close();
		# consensus_stats.PlotConsensusStats(dataset_name, consensus_prefix, consensus_statistics, out_consensus_plot_png);

	if (mode_code & MODE_CODE_CALC_ERROR_RATE):
		sys.stderr.write('[analyzesam.py] Calculating error rates - individual indels...\n');
		errorrates.ProcessFromFiles(reference_file, sam_file, out_accuracy_counts_path, False);
		sys.stderr.write('\n');
	if ((mode_code & MODE_CODE_CALC_ERROR_RATE) or (mode_code & MODE_CODE_COLLECT_ERROR_RATE)):
		# Collecting stats from the error rate estimation.
		error_rates_return = errorrates.CollectAccuracy(sam_file, out_accuracy_counts_path, False);

		try:
			[summary_lines, summary_cigar, error_rate_hist, insertion_hist, deletion_hist, snp_hist, match_hist, cigar_dataset_name, out_png_path] = error_rates_return;
			fp.write('[CIGAR statistics - individual indels]\n');
			fp.write(summary_lines + '\n');
			sys.stderr.write('[CIGAR statistics - individual indels]\n');
			sys.stderr.write(summary_lines + '\n');
			summary += '[CIGAR statistics - individual indels]\n';
			summary += summary_lines + '\n\n';

		except Exception, e:
			sys.stderr.write(str(e) + '\n');
			sys.stderr.write('Returned values: %s\n' % (str(error_rates_return)));

	if (mode_code & MODE_CODE_CALC_ERROR_RATE_INDELS_AS_EVENTS):
		sys.stderr.write('[analyzesam.py] Calculating error rates - indels as events...\n');
		errorrates.ProcessFromFiles(reference_file, sam_file, out_accuracy_counts_indel_events_path, True);
		sys.stderr.write('\n');
	if ((mode_code & MODE_CODE_CALC_ERROR_RATE_INDELS_AS_EVENTS) or (mode_code & MODE_CODE_COLLECT_ERROR_RATE_INDELS_AS_EVENTS)):
		# Collecting stats from the error rate estimation.
		error_rates_return = errorrates.CollectAccuracy(sam_file, out_accuracy_counts_indel_events_path, False);
		
		try:
			[summary_lines, summary_cigar, error_rate_hist, insertion_hist, deletion_hist, snp_hist, match_hist, cigar_dataset_name, out_png_path] = error_rates_return;
			# [summary_cigar, error_rate_hist, insertion_hist, deletion_hist, snp_hist, match_hist, cigar_dataset_name, out_png_path] = errorrates.CollectAccuracy(sam_file, out_accuracy_counts_indel_events_path, False);
			fp.write('[CIGAR statistics - indels as events]\n');
			fp.write(summary_lines + '\n');
			sys.stderr.write('[CIGAR statistics - indels as events]\n');
			sys.stderr.write(summary_lines + '\n');
			summary += '[CIGAR statistics - indels as events]\n';
			summary += summary_lines + '\n\n';
		except Exception, e:
			sys.stderr.write(str(e) + '\n');
			sys.stderr.write('Returned values: %s\n' % (str(error_rates_return)));

	# if (mode_code & MODE_CODE_CALC_READ_COUNTS):
	# 	sys.stderr.write('[main_sam_analysis.py] Counting mapped reads...\n');
	# 	[num_alignments, num_mapped_alignments, num_unique_reads, num_mapped_reads, num_mapped_bases] = utility_sam.CountMappedReads(sam_file);
	# 	sys.stderr.write('\n');
	# 	summary +=

	num_alignments = 0;
	num_mapped_alignments = 0;
	num_unique_reads = 0;
	num_mapped_reads = 0;
	num_mapped_bases = 0;
	if ((mode_code & MODE_CODE_CALC_READ_COUNTS) or (mode_code & MODE_CODE_COLLECT_READ_COUNTS)):
		sys.stderr.write('[main_sam_analysis.py] Counting mapped reads...\n');
		[num_alignments, num_mapped_alignments, num_unique_reads, num_mapped_reads, num_mapped_bases] = utility_sam.CountMappedReads(sam_file);
		sys.stderr.write('\n');
		summary_read_count = '';
		summary_read_count += 'num_alignments: %d\n' % num_alignments;
		summary_read_count += 'num_mapped_alignments: %d (%.2f%%)\n' % (num_mapped_alignments, (float(num_mapped_alignments) / float(num_alignments)) * 100.0);
		summary_read_count += 'num_unmapped_alignments: %d (%.2f%%)\n' % ((num_alignments - num_mapped_alignments), (float(num_alignments - num_mapped_alignments) / float(num_alignments)) * 100.0);
		summary_read_count += 'num_mapped_reads: %d\n' % num_mapped_reads;
		summary_read_count += 'num_uniquely_mapped_reads: %d\n' % num_unique_reads;
		summary_read_count += 'num_mapped_bases: %d\n' % num_mapped_bases;
		summary_read_count += 'num_read_in_input_reads_file: %d\n' % (fastqinfo_num_seqs);
		summary_read_count += 'num_bases_in_input_reads_file: %d\n' % (fastqinfo_total_seq_len);
		summary_read_count += 'percent_mapped_reads: %.2f%%\n' % ((float(num_mapped_reads) / float(fastqinfo_num_seqs)) * 100.0);
		summary_read_count += 'percent_mapped_bases: %.2f%%\n' % ((float(num_mapped_bases) / float(fastqinfo_total_seq_len)) * 100.0);

		summary += '[SAM file statistics]\n';
		summary += summary_read_count + '\n\n';

		###############################3
		# This version is more verbose, but is different than figure1.xlsx!
		# csv_line += '%d\t%d\t%f\t%f\t' % (num_mapped_reads, num_mapped_bases, (float(num_mapped_reads)/float(fastqinfo_num_seqs))*100.0, (float(num_mapped_bases)/float(fastqinfo_total_seq_len))*100.0);
		# csv_header += 'num_mapped_reads\tnum_mapped_bases\tpercent_mapped_reads\tpercent_mapped_bases\t';
		###############################3

		###############################3
		# This version is in the same order as figure1.xlsx!
		csv_line += '%d\t%f\t%d\t%f\t' % (num_mapped_bases, (float(num_mapped_bases)/float(fastqinfo_total_seq_len))*100.0, num_mapped_reads, (float(num_mapped_reads)/float(fastqinfo_num_seqs))*100.0);
		csv_header += 'Num mapped bases\tPercent mapped bases\tNum mapped reads\tPercent mapped reads\t';
		###############################3

		# Collecting stats from the SAM file.
		# [summary_read_count_readable, read_count_plot_lines] = count_mapped_reads.CollectSAMStats(out_count_mapped_reads_prefix, suppress_error_messages);
		sys.stderr.write('[SAM file statistics]\n');
		sys.stderr.write(summary_read_count + '\n');
		fp.write('[SAM file statistics]\n');
		fp.write(summary_read_count + '\n');

	# Get the CPU time and memory consumption from the process.
	[cmdline, realtime, cputime, usertime, systemtime, maxrss, rsscache, memtime_lines] = ParseMemTime(sam_file);
	reads_per_sec = 0.0 if (cputime == 0) else (float(num_unique_reads)/float(cputime));
	bases_per_sec = 0.0 if (cputime == 0) else (float(num_mapped_bases)/float(cputime));
	bases_per_mb = 0.0 if (cputime == 0) else float(num_mapped_bases)/float(maxrss);

	summary_memtime = '\n'.join(memtime_lines);
	summary += '[Memtime]\n' + summary_memtime + '\n\n';
	csv_line += '%f\t%f\t' % (cputime, maxrss);
	csv_line += '%f\t%f\t%f\t' % (reads_per_sec, bases_per_sec, bases_per_mb);
	csv_header += 'CPU time [s]\tMemory [MB]\treads/sec\tbases/sec\tbases/MB\t';
	# Just simply verbose to screen and summary file.
	sys.stderr.write('[Memtime]\n');
	sys.stderr.write(summary_memtime + '\n\n');
	fp.write('[Memtime]\n');
	fp.write(summary_memtime + '\n\n');

	# sys.stderr.write('[SAM headers]\n');
	# sys.stderr.write(summary_file_sam_headers + '\n');
	# sys.stderr.write('[SAM file statistics]:\n');
	# sys.stderr.write(summary_read_count_readable + '\n');
	# sys.stderr.write('[CIGAR statistics]\n');
	# sys.stderr.write(summary_cigar + '\n');
	# sys.stderr.write('[Consensus statistics]\n');
	# sys.stderr.write(summary_consensus + '\n');
	# sys.stderr.write('\n');
	sys.stderr.write(summary);
	sys.stderr.write('[Done!]\n');
	sys.stderr.write('\n');

	fp.close();

	csv_header += 'Coverage threshold';
	csv_line += '%d\t' % (consensus_coverage_threshold);

	csv_header = csv_header.strip();
	csv_line = csv_line.strip();

	if (mode_code & MODE_CODE_HEADER):
		return '%s\n%s\n' % (csv_header, csv_line);

	return '%s\n' % (csv_line);
	
	#plt.show();

def ReadlinesWrapper(file_path):
	try:
		fp = open(file_path, 'r');
		lines = fp.readlines();
		fp.close();
		return [line.strip() for line in lines];
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!' % file_path);
	
	return [];

def WritelineWrapper(line, file_path):
	try:
		fp = open(file_path, 'w');
		fp.write(line);
		fp.close();
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!' % file_path);

# def ReformatLinesConsensusForSingleCoverage(lines_consensus, coverage):
# 	if (len(lines_consensus) < 3):
# 		print 'ERROR: Not enough lines in the input file!';
# 		return lines_consensus;
	
# 	split_line = lines_consensus[2].split('\t');
# 	if (coverage in split_line):
# 		column_index = split_line.index(coverage);
# 	else:
# 		print 'ERROR: Coverage "%s" is not in the input file!' % str(coverage);
# 		return [];
	
# 	ret_lines = [];
# 	mapper_dict = {};
	
# 	i = 3;
# 	while (i < len(lines_consensus)):
# 		split_line = lines_consensus[i].split('\t');
# 		mapper_name = split_line[1];
# 		param_name = split_line[2];
# 		param_value = split_line[column_index];
		
# 		try:
# 			param_dict = mapper_dict[mapper_name];
# 			param_dict[param_name] = param_value;
# 		except:
# 			mapper_dict[mapper_name] = {param_name: param_value};
		
# 		i += 1;
	
# 	ret_lines.append(lines_consensus[0]);
# 	ret_lines.append(lines_consensus[1]);
	
# 	final_param_names = None;
# 	for mapper_name in mapper_dict.keys():
# 		#print mapper_name, mapper_dict[mapper_name];
# 		param_dict = mapper_dict[mapper_name];
# 		param_names = sorted(param_dict.keys());
# 		param_values = [param_dict[param_name] for param_name in param_names];
		
# 		if (final_param_names == None):
# 			final_param_names = param_names;
# 			ret_lines.append('mapper_name\t' + '\t'.join(param_names));
# 		else:
# 			if (param_names != final_param_names):
# 				print 'ERROR: Parameters are not the same for different mappers!';
# 				return [];
		
# 		ret_lines.append(mapper_name + '\t' + '\t'.join(param_values));
	
# 	return ret_lines;
	
# def CollectAndFormatAll(folder_path, sam_files, suppress_error_messages=False, consensus_single_coverage=''):
# 	# Create the output path if it doesn't exist yet.
# 	output_folder = folder_path + '/analysis-final';
# 	if not os.path.exists(output_folder):
# 		os.makedirs(output_folder);
	
# 	all_accuracies = [];
# 	all_consensus = [];
# 	all_consensus_single_coverage = [];
# 	all_read_counts = [];
# 	all_memtimes = [];

# 	intermediate_output_folder = folder_path + '/analysis-intermediate';
	
# 	for sam_file in sam_files:
# 		sam_basename = os.path.splitext(os.path.basename(sam_file))[0];
# 		# Path to intermediate results.
# 		intermediate_accuracy_counts_plot_path = '%s/accuracy-%s.csv.plot' % (intermediate_output_folder, sam_basename);
# 		intermediate_consensus_prefix = '%s/consensus-%s' % (intermediate_output_folder, sam_basename);
# 		intermediate_consensus_plot_path = intermediate_consensus_prefix + '.plot';
# 		intermediate_count_mapped_reads_plot_path = '%s/count_reads-%s.plot' % (intermediate_output_folder, sam_basename);
		
# 		lines_accuracy = ReadlinesWrapper(intermediate_accuracy_counts_plot_path);
# 		lines_consensus = ReadlinesWrapper(intermediate_consensus_plot_path);
# 		lines_read_counts = ReadlinesWrapper(intermediate_count_mapped_reads_plot_path);
		
# 		lines_consensus_single_coverage = [];
# 		if (consensus_single_coverage != ''):
# 			lines_consensus_single_coverage = ReformatLinesConsensusForSingleCoverage(lines_consensus, consensus_single_coverage);
		
# 		if (len(lines_accuracy) > 3):
# 			all_accuracies += (lines_accuracy[2:] if (len(all_accuracies) == 0) else lines_accuracy[3:]);
# 		if (len(lines_consensus) > 3):
# 			all_consensus += (lines_consensus[2:] if (len(all_consensus) == 0) else lines_consensus[3:]);
# 		if (len(lines_consensus_single_coverage) > 3):
# 			all_consensus_single_coverage += (lines_consensus_single_coverage[2:] if (len(all_consensus_single_coverage) == 0) else lines_consensus_single_coverage[3:]);
# 		if (len(lines_read_counts) > 3):
# 			all_read_counts += (lines_read_counts[2:] if (len(all_read_counts) == 0) else lines_read_counts[3:]);
		
# 		[lines_memtime, lines_formatted_memtime] = memtime_stats.CollectMemTime(sam_file);
# 		if (len(lines_formatted_memtime) > 1):
# 			all_memtimes += (lines_formatted_memtime[0:] if (len(all_memtimes) == 0) else lines_formatted_memtime[1:]);
	
# 	dataset_name = os.path.basename(folder_path);
# 	out_final_accuracy = '%s/plot-accuracy-%s.csv' % (output_folder, dataset_name);
# 	out_final_consensus = '%s/plot-consensus-%s.csv' % (output_folder, dataset_name);
# 	if (consensus_single_coverage != ''):
# 		out_final_consensus_single_coverage = '%s/plot-consensus-%s-cov_%s.csv' % (output_folder, dataset_name, consensus_single_coverage);
# 	out_final_read_counts = '%s/plot-read_counts-%s.csv' % (output_folder, dataset_name);
# 	out_final_memtimes = '%s/plot-memtimes-%s.csv' % (output_folder, dataset_name);
	
# 	WritelineWrapper('\n'.join(all_accuracies), out_final_accuracy);
# 	WritelineWrapper('\n'.join(all_consensus), out_final_consensus);
# 	WritelineWrapper('\n'.join(all_consensus_single_coverage), out_final_consensus_single_coverage);
# 	WritelineWrapper('\n'.join(all_read_counts), out_final_read_counts);
# 	WritelineWrapper('\n'.join(all_memtimes), out_final_memtimes);



# def ProcessSamFiles(sam_files, reference_file, output_folder_path, num_reads_in_reference_sam=1, suppress_error_messages=False, filter_filename=''):
# 	i = 0;
# 	for sam_file in sam_files:
# 		i += 1;
		
# 		if (len(filter_filename) > 0 and (not (filter_filename in sam_file))):
# 			continue;
		
# 		Main(sam_file, reference_file, output_folder_path, num_reads_in_reference_sam, suppress_error_messages, [20]);

# 	CollectAndFormatAll(output_folder_path, sam_files, suppress_error_messages, '20');

# def MainFromFolder(folder_path, reference_file, num_reads_in_reference_sam=1, suppress_error_messages=False, filter_filename=''):
# 	sam_files = filesandfolders.FindFiles(folder_path, '*.sam');
	
# 	ProcessSamFiles(sam_files, reference_file, folder_path, num_reads_in_reference_sam, suppress_error_messages, filter_filename);




# Example usage:
#	./main_sam_analysis.py file ../data/minion-review/reads-E.Coli-R7.3/alignment/ecoliR7.3/graphmap-params_experimenting-newest-filtered.sam ../data/minion-review/reference/escherichia_coli.fa
#	./main_sam_analysis.py file ../reads-simulated/OxfordNanopore-pbsim-30_percent-100k/saccharomyces_cerevisiae/reads-1000.sam ../reference-genomes/saccharomyces_cerevisiae.fa

# Example usage:
#	./analyzesam.py [file|folder] <sam_path> <reference_path> <reads_path> [<simulated_sam_path>]
if __name__ == "__main__":
	if (len(sys.argv) >= 6 and len(sys.argv) <= 8):
		run_from_args(sys.argv[1:]);
	else:
		print_usage_and_exit();
