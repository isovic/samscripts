#! /usr/bin/env python

# Copyright Ivan Sovic, 2015. www.sovic.org
#
# Creates a pileup from a given SAM/BAM file, and calls consensus of structural variants from the mappings.

import os;
import sys;
import operator;
import subprocess;

class StructuralEvent:
	def __init__(self):
		self.clear();
		
	def clear(self):
		self.type = '';
		self.start = -1;
		self.end = -1;
		self.open_coverage = -1;
		self.close_coverage = -1;
		self.max_coverage = -1;
		self.max_event_length = -1;
		self.comment = 'real';

	def is_open(self):
		if (self.start >= 0 and self.end == -1):
			return True;
		return False;

	def to_string(self):
		# ret = 'type: %s\tstart: %d\tend: %d\topen_cov: %d\tclose_cov: %d\tmax_cov: %d' % (self.type, self.start, self.end, self.open_coverage, self.close_coverage, self.max_coverage);
		ret = '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s' % (self.type, self.start, self.end, (self.end - self.start), self.open_coverage, self.close_coverage, self.max_coverage, self.max_event_length, self.comment);
		return ret;

	def csv_header(self):
		ret = 'type\tstart\tend\tlength\topen_cov\tclose_cov\tmax_cov\tmax_event_length\tcomment';
		return ret;

def increase_in_dict(dict_counter, value):
	try:
		dict_counter[value] += 1;
	except:
		dict_counter[value] = 1;

def process_mpileup_line(line, line_number, ret_variant_list, ret_vcf_list, ret_snp_count, ret_insertion_count, ret_deletion_count, ret_num_undercovered_bases, ret_num_called_bases, ret_num_correct_bases, ret_coverage_sum, ret_opened_insertion_events, ret_opened_deletion_events, coverage_threshold, verbose=False):
	# Split the line, and perform a sanity check.
	split_line = line.strip().split('\t');
	if (len(split_line) < 5 or len(split_line) > 6):
		sys.stderr.write(line + '\n');
		return [0, 0];

	ref_name = split_line[0];
	position = split_line[1];
	ref_base = split_line[2];
	coverage = split_line[3];
	original_bases = split_line[4];
	if (len(split_line) == 6):
		qualities = split_line[5];
	
	bases = '';
	
	# Replace the '.' and ',' signs with the actual reference base.
	i = 0;
	while (i < len(original_bases)):
		if (original_bases[i] == '.' or original_bases[i] == ','):
			bases += ref_base;
		else:
			bases += original_bases[i];
		i += 1;

	base_counts = {};
	insertion_count = 0;
	current_base_deletion_count = 0;
	deletion_count = 0;
	insertion_event_counts = {};
	deletion_event_counts = {};
	end_counts = 0;

	# print 'position: %s' % position;
	# print 'bases: "%s"' % bases;
	# print 'line_number: %d' % line_number;
	# print line;
	# print '';
	# sys.stdout.flush();

	i = 0;
	while (i < len(bases)):
		base = bases[i];
		
		if (base == r'^'):
			# This is the starting position of a read. It encodes two
			# symbols: '^' marking the read start and a char marking the
			# mapping quality of the read.
			#increase_in_dict(base_counts, bases[i + 1].upper());
			i += 1;			# Increase only by 1, because we have i += 1 down there.
		elif (base == r'$'):
			# This marks the end of a read.
			end_counts += 1;
		elif (base == r'*'):
			# This is a deletion, just count it.
			current_base_deletion_count += 1;
		elif (base == r'-'):
			# This marks the occurance of deletions. It is a composite object
			# consisting of: the special character '-', the number of the deleted bases
			# and the actual bases that are deleted (these bases follow the current position).
			# In our approach, we ignore this case, because we count deletions one by one
			# through the '*' character.
			
			# Get the number of bases that need to be skipped in the string.
			j = (i + 1);
			while (bases[j] in '0123456789'):
				j += 1;
			num_bases = int(bases[(i + 1):j]);
			skip_bases = (j - i) + num_bases - 1;
			deletion_count += 1;
			deletion = bases[j : (j + num_bases)].upper();
			increase_in_dict(deletion_event_counts, deletion);
			# Skip the length of the numeric entry plus the actual number of bases
			# that need to be skipped.
			i += skip_bases;
		elif (base == r'+'):
			# This marks the occurance of an insertion. It is a composite object
			# consisting of: the special character '+', the number of the inserted bases
			# and the actual bases that are inserted (these bases follow the current position).
			# Similar to the deletion marking, but here we actually care about the bases,
			# and we need to make an allele aware count.

			# Get the number of bases that are inserted;
			j = (i + 1);
			while (bases[j] in '0123456789'):
				j += 1;
			num_bases = int(bases[(i + 1):j]);
			skip_bases = (j - i) + num_bases - 1;
			insertion_count += 1;
			insertion = bases[j : (j + num_bases)].upper();
			increase_in_dict(insertion_event_counts, insertion);
			i += skip_bases;
		else:
			increase_in_dict(base_counts, bases[i].upper());
		i += 1;
	
	event_length_threshold = 100;
	# event_length_threshold = 50;
	# event_length_threshold = 10;
	event_length_threshold = 20;
	# event_length_threshold = 9;
	# event_length_threshold = 5;
	# event_coverage_threshold = 10;

	# sorted_insertion_counts = sorted(insertion_event_counts.items(), key=operator.itemgetter(1));
	for insertion_event in insertion_event_counts.keys():
		if (len(insertion_event) > event_length_threshold):
			# print '';
			# print '\n\tInsertion event length: %d, position: %d, len(opened_insertion_events): %d' % (len(insertion_event), int(position), len(ret_opened_insertion_events));
			ret_opened_insertion_events.append([len(insertion_event), len(insertion_event)]);

	for deletion_event in deletion_event_counts.keys():
		if (len(deletion_event) > event_length_threshold):
			# print '\n\tDeletion event length: %d, position: %d, len(opened_deletion_events): %d' % (len(deletion_event), int(position), len(ret_opened_deletion_events));
			ret_opened_deletion_events.append([len(deletion_event), len(deletion_event)]);

	return [position, coverage];

def process_mpileup(alignments_path, reference_path, mpileup_path, coverage_threshold, out_file, thread_id=0, bed_position=''):
	fp = None;
	try:
		fp = open(mpileup_path, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % mpileup_path);
		return;
	
	ret_variant_list = [];
	ret_vcf_list = [];
	ret_snp_count = [0];
	ret_insertion_count = [0];
	ret_deletion_count = [0];
	ret_num_undercovered_bases = [0];
	ret_num_called_bases = [0];
	ret_num_correct_bases = [0];
	ret_coverage_sum = [0];
	
	lines = fp.readlines();
	fp.close();
	
	fp_variant = None;
	fp_vcf = None;
	if (not os.path.exists(os.path.dirname(out_file))):
		sys.stderr.write('Creating the output folder: "%s".\n' % (os.path.dirname(out_file)));
		os.makedirs(os.path.dirname(out_file));

	fp_out = open(out_file, 'w');
	# if (output_prefix != ''):
		# variant_file = ('%s-cov_%d.variant.csv' % (output_prefix, coverage_threshold));
		# fp_variant = open(variant_file, 'w');

		# vcf_file = ('%s-cov_%d.variant.vcf' % (output_prefix, coverage_threshold));
		# fp_vcf = open(vcf_file, 'w');
		# fp_vcf.write('##fileformat=VCFv4.0\n');
		# fp_vcf.write('##fileDate=20150409\n');
		# fp_vcf.write('##source=%s\n' % (' '.join(sys.argv)));
		# fp_vcf.write('##reference=%s\n' % reference_path);
		# fp_vcf.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">\n');
		# fp_vcf.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n');
		# fp_vcf.write('##INFO=<ID=SB,Number=1,Type=Integer,Description="Phred-scaled strand bias at this position">\n');
		# fp_vcf.write('##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">\n');
		# fp_vcf.write('##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">\n');
		# fp_vcf.write('##INFO=<ID=CONSVAR,Number=0,Type=Flag,Description="Indicates that the variant is a consensus variant (as opposed to a low frequency variant).">\n');
		# fp_vcf.write('##INFO=<ID=HRUN,Number=1,Type=Integer,Description="Homopolymer length to the right of report indel position">\n');
		# fp_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n');
		# fp_vcf.flush();



	use_bed = False;
	bed_chromosome = "";
	bed_pos_start = 0;
	bed_pos_end = len(lines);
	if (bed_position != ""):
		bed_split = bed_position.split(':');
		if (len(bed_split) != 2):
			use_bed = False;
		else:
			bed_chromosome = bed_split[0];
			bed_pos_split = bed_split[1].split('-');
			if (len(bed_pos_split) != 2):
				use_bed = False;
			else:
				bed_pos_start = int(bed_pos_split[0]);
				bed_pos_end = int(bed_pos_split[1]);
				use_bed = True;
				sys.stderr.write('Using location specified through commandline:\n');
				sys.stderr.write('\tChromosome: "%s"\n' % bed_chromosome);
				sys.stderr.write('\tStart: %d\n' % bed_pos_start);
				sys.stderr.write('\tEnd: %d\n\n' % bed_pos_end);

	opened_insertion_events = [];
	opened_deletion_events = [];
	threshold_switch_insertions = False;
	threshold_switch_deletions = False;
	num_insertion_events = 0;
	num_deletion_events = 0;

	mpileup_pos = 0;
	mpileup_cov = 0;

	new_insertion_event = StructuralEvent();
	new_deletion_event = StructuralEvent();
	insertion_events = [];
	deletion_events = [];

	fp_out.write(new_insertion_event.csv_header() + '\n');

	# i = 0;
	i = 0 if (use_bed == False) else max((bed_pos_start - 10), 0);
	j = 0;
	while (i < bed_pos_end): # len(lines)):
		line = lines[i];

		if (use_bed == True):
			line_split = line.strip().split('\t');
			if (len(line_split) > 2 and line_split[0] == bed_chromosome):
				current_pos = int(line_split[1]);
				if (current_pos < bed_pos_start or current_pos >= bed_pos_end):
					i += 1;
					j += 1;
					continue;
			else:
				# print line_split[0];
				# print bed_chromosome;
				i += 1;
				j += 1;
				continue;

		if (thread_id == 0):
			if ((j % 1000) == 0):
				sys.stderr.write('\r[%d] insertion events: %d, deletion events: %d, current coverage: %d' % (i, num_insertion_events, num_deletion_events, int(mpileup_cov)));
				sys.stderr.flush();

		variant_list_length = len(ret_variant_list);
		vcf_list_length = len(ret_vcf_list);

		[mpileup_pos, mpileup_cov] = process_mpileup_line(line, i, ret_variant_list, ret_vcf_list, ret_snp_count, ret_insertion_count, ret_deletion_count,
														  ret_num_undercovered_bases, ret_num_called_bases, ret_num_correct_bases, ret_coverage_sum,
														  opened_insertion_events, opened_deletion_events, coverage_threshold, verbose=use_bed);

		event_coverage_threshold = max(int(0.15 * float(mpileup_cov)), 5);

		if (threshold_switch_insertions == False and len(opened_insertion_events) >= event_coverage_threshold):
			##### try:
			##### 	sys.stderr.write('\n+++ [I] Opened plausible structural event! Position: %d, coverage: %d, len(opened_insertion_events): %d, last_len: %d, event_coverage_threshold: %d\n' % (int(mpileup_pos), int(mpileup_cov), len(opened_insertion_events), opened_insertion_events[-1][0], event_coverage_threshold));
			##### except Exception, e:
			##### 	sys.stderr.write(str(e));

			num_insertion_events += 1;
			threshold_switch_insertions = True;

			new_insertion_event.start = int(mpileup_pos);
			new_insertion_event.type = 'insertion_in_read';
			new_insertion_event.open_coverage = len(opened_insertion_events);
			new_insertion_event.close_coverage = -1;
			new_insertion_event.max_coverage = len(opened_insertion_events);
			new_insertion_event.max_event_length = max([new_insertion_event.max_event_length] + [insertion_event[1] for insertion_event in opened_insertion_events]);
			new_insertion_event.end = -1;
		elif (threshold_switch_insertions == True and len(opened_insertion_events) < event_coverage_threshold):
			##### sys.stderr.write('\n--- [I] Closed plausible structural event! Position: %d, coverage: %d, len(opened_insertion_events): %d\n' % (int(mpileup_pos), int(mpileup_cov), len(opened_insertion_events)));
			threshold_switch_insertions = False;
			new_insertion_event.end = int(mpileup_pos)
			new_insertion_event.close_coverage = len(opened_insertion_events);
			insertion_events.append(new_insertion_event);
			##### sys.stderr.write('Event:\n');
			sys.stderr.write(' ' + new_insertion_event.to_string());
			sys.stderr.write('\n');
			fp_out.write(new_insertion_event.to_string() + '\n');
			fp_out.flush();
			new_insertion_event.clear();

		if (threshold_switch_deletions == False and len(opened_deletion_events) >= event_coverage_threshold):
			##### try:
			##### 	sys.stderr.write('\n+++ [D] Opened plausible structural event! Position: %d, coverage: %d, len(opened_deletion_events): %d, last_len: %d, event_coverage_threshold: %d\n' % (int(mpileup_pos), int(mpileup_cov), len(opened_deletion_events), opened_deletion_events[-1][0], event_coverage_threshold));
			##### except Exception, e:
			##### 	sys.stderr.write(str(e));
			threshold_switch_deletions = True;
			num_deletion_events += 1;
			
			new_deletion_event.start = int(mpileup_pos);
			new_deletion_event.type = 'deletion_in_read';
			new_deletion_event.open_coverage = len(opened_deletion_events);
			new_deletion_event.close_coverage = -1;
			new_deletion_event.max_coverage = len(opened_deletion_events);
			new_deletion_event.max_event_length = max([new_deletion_event.max_event_length] + [deletion_event[1] for deletion_event in opened_deletion_events]);
			new_deletion_event.end = -1;
		elif (threshold_switch_deletions == True and len(opened_deletion_events) < event_coverage_threshold):
			##### sys.stderr.write('\n--- [D] Closed plausible structural event! Position: %d, coverage: %d, len(opened_deletion_events): %d\n' % (int(mpileup_pos), len(mpileup_cov), len(opened_deletion_events)));
			threshold_switch_deletions = False;
			new_deletion_event.end = int(mpileup_pos)
			new_deletion_event.close_coverage = len(opened_insertion_events);
			deletion_events.append(new_deletion_event);
			##### sys.stderr.write('Event:\n');
			sys.stderr.write(' ' + new_deletion_event.to_string());
			sys.stderr.write('\n');
			fp_out.write(new_deletion_event.to_string() + '\n');
			fp_out.flush();
			new_deletion_event.clear();

		if (new_insertion_event.is_open() == True):
			new_insertion_event.max_coverage = max(new_insertion_event.max_coverage, len(opened_insertion_events));
			new_insertion_event.max_event_length = max([new_insertion_event.max_event_length] + [insertion_event[1] for insertion_event in opened_insertion_events]);
		if (new_deletion_event.is_open() == True):
			new_deletion_event.max_coverage = max(new_deletion_event.max_coverage, len(opened_deletion_events));
			new_deletion_event.max_event_length = max([new_deletion_event.max_event_length] + [deletion_event[1] for deletion_event in opened_deletion_events]);

		# event_coverage_threshold = 10;

		# if (len(opened_insertion_events) > coverage_threshold):
		# 	print 'Plausible structural event! Position: %d, coverage: %d' % (i, len(opened_insertion_events));

		# if (len(opened_insertion_events) > 0):
		# 	print opened_insertion_events;
		i1 = 0;
		while (i1 < len(opened_insertion_events)):
			opened_insertion_events[i1][0] -= 1;
			i1 += 1;
		opened_insertion_events = filter(lambda a: a[0] > 0, opened_insertion_events);
		# if (len(opened_insertion_events) > 0):
		# 	print opened_insertion_events;
		# 	print '';

		i1 = 0;
		while (i1 < len(opened_deletion_events)):
			opened_deletion_events[i1][0] -= 1;
			i1 += 1;
		opened_deletion_events = filter(lambda a: a[0] > 0, opened_deletion_events);



		
		# if (len(ret_variant_list) > variant_list_length and fp_variant != None):
		# 	fp_variant.write('\n'.join(ret_variant_list[variant_list_length:]) + '\n');
		# 	fp_variant.flush();

		# if (len(ret_vcf_list) > vcf_list_length and fp_vcf != None):
		# 	fp_vcf.write('\n'.join(ret_vcf_list[vcf_list_length:]) + '\n');
		# 	fp_vcf.flush();

		# i += num_bases_to_skip;
		i += 1;
		j += 1;
		
		#if (i > 10000):
			#break;

	# if (fp_variant != None):
	# 	fp_variant.close();

	# if (fp_vcf != None):
	# 	fp_vcf.close();
	
	sys.stderr.write('\n');

	# summary_lines = '';
	# summary_lines += 'alignments_file: %s\n' % alignments_path;
	# summary_lines += 'mpileup_file: %s\n' % mpileup_path;
	# summary_lines += 'coverage_threshold: %d\n' % coverage_threshold;
	# summary_lines += 'snp_count: %d\n' % ret_snp_count[0];
	# summary_lines += 'insertion_count: %d\n' % ret_insertion_count[0];
	# summary_lines += 'deletion_count: %d\n' % ret_deletion_count[0];
	# summary_lines += 'num_undercovered_bases: %d\n' % ret_num_undercovered_bases[0];
	# summary_lines += 'num_called_bases: %d\n' % ret_num_called_bases[0];
	# summary_lines += 'num_correct_bases: %d\n' % ret_num_correct_bases[0];
	# summary_lines += 'average_coverage: %.2f\n' % ((float(ret_coverage_sum[0])/float((i + 1))));
	
	# sys.stderr.write(summary_lines + '\n');
	sys.stderr.write('\n');
	fp_out.close();
	
	# if (output_prefix != ''):
	# 	#summary_file = output_prefix + '.conssum';
	# 	summary_file = ('%s-cov_%d.variant.sum' % (output_prefix, coverage_threshold));

	# 	try:
	# 		fp_sum = open(summary_file, 'w');
	# 		fp_sum.write(summary_lines);
	# 		fp_sum.close();
	# 	except IOError:
	# 		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % (summary_file));

class StructEvent:
	def __init__(self, size=0, type='', position=0):
		self.size = size;
		self.type = type;
		self.position = position;
		self.position_final = position;
		self.open_cov = 0;
		self.close_cov = 0;
		self.max_cov = 0;
		self.max_event_length = 0;
		self.comment = '';
		self.evaluation = '(not evaluated)';

	def verbose(self):
		sys.stderr.write(self.verbose_to_string());

	def verbose_to_string(self):
		event_type = 'insertion_in_read' if (self.type == 'D' or self.type == 'insertion_in_read') else 'deletion_in_read';
		return 'Event type: %s, position: %7d, end: %7d, size: %4d, max_event_length = %4d, result: %s, comment: %s' % (event_type, self.position, (self.position + self.size), self.size, self.max_event_length, self.evaluation, self.comment);

	def csv_line(self):
		# ret = 'type: %s\tstart: %d\tend: %d\topen_cov: %d\tclose_cov: %d\tmax_cov: %d' % (self.type, self.start, self.end, self.open_coverage, self.close_coverage, self.max_coverage);
		event_type = 'insertion_in_read' if (self.type == 'D' or self.type == 'insertion_in_read') else 'deletion_in_read';
		ret = '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s' % (event_type, self.position, (self.position + self.size), self.size, self.open_cov, self.close_cov, self.max_cov, self.max_event_length, self.comment);
		return ret;

	def csv_header(self):
		ret = 'type\tstart\tend\tlength\topen_cov\tclose_cov\tmax_cov\tmax_event_length\tcomment';
		return ret;

def chain_events(events_path):
	events = [];

	num_lines = 0;
	fp = open(events_path, 'r');
	for line in fp:
		num_lines += 1;
		line = line.strip();
		if (len(line) == 0):
			continue;

		if (num_lines == 1):
			labels = line;
			continue;

		# print line;

		values = line.split('\t');
		# event_type = 'I' if (values[0] == 'deletion_in_read') else 'D';
		event_type = 'insertion_in_read' if (values[0] == 'D' or values[0] == 'insertion_in_read') else 'deletion_in_read';
		start = int(values[1]);
		end = int(values[2]);
		length = int(values[3]);
		open_cov = int(values[4]);
		close_cov = int(values[5]);
		max_cov = int(values[6]);
		max_event_length = int(values[7]);
		comment = values[8] if (len(values) > 8) else '';

		event = StructEvent(length, event_type, start);
		event.open_cov = open_cov;
		event.close_cov = close_cov;
		event.max_cov = max_cov;
		event.max_event_length = max_event_length;
		event.comment = comment;

		events.append(event);

	fp.close();

	### Chain events that got fragmented
	sorted_events = sorted(events, key=lambda x: x.position);
	chained_events = [];
	for event in sorted_events:
		if (len(chained_events) > 0):
			last_size = chained_events[-1].size;
			last_endpoint = chained_events[-1].position + last_size;
			last_type = chained_events[-1].type;
			distance = abs(event.position - last_endpoint);
			if (distance <= (last_size + event.size) and last_type == event.type):
				chained_events[-1].close_cov = event.close_cov;
				chained_events[-1].max_cov = max(chained_events[-1].max_cov, event.max_cov);
				chained_events[-1].max_event_length = max(chained_events[-1].max_event_length, event.max_event_length);
				chained_events[-1].size = event.position + event.size - chained_events[-1].position + 1;
				# print 'Chaining: %s' % chained_events[-1].verbose_to_string();
			else:
				chained_events.append(event);
		else:
			chained_events.append(event);

	# sys.stderr.write('Chained events:\n');
	# for event in chained_events:
	# 	sys.stderr.write(event.csv_line() + '\n');

	fp = open(events_path + '.chained.csv', 'w');
	fp.write(labels + '\n');
	fp.write('\n'.join([event.csv_line() for event in chained_events]));
	fp.write('\n');
	fp.close();

# Problem je sto je zadnji entry u chainanim eventovima falio u output! Treba provjeriti gdje se proguta!

	return chained_events;
	# return events;

def main(alignments_path, reference_path, coverage_threshold, out_file, thread_id=0, bed_position=""):
	# Sanity checking the existence of the file, and the correctness of its extension.
	# Also, if input file is a SAM file, then convert it to a sorted BAM.
	alignments_path_bam = alignments_path;
	if (os.path.exists(alignments_path) == False):
		sys.stderr.write('ERROR: File "%s" does not exist!\n' % alignments_path);
		return;
	if (alignments_path.endswith('sam')):
		# Determine the path where the new BAM file will be generated.
		dir_name = os.path.dirname(alignments_path);
		if (dir_name == ''):
			dir_name = '.';
		alignments_path_bam = dir_name + '/' + os.path.splitext(os.path.basename(alignments_path))[0] + '.bam'
		alignments_path_bam_exists = os.path.exists(alignments_path_bam);
		# Check if a BAM file with the given name already exists.
		if (alignments_path_bam_exists == False or (alignments_path_bam_exists == True and os.path.getmtime(alignments_path) > os.path.getmtime(alignments_path_bam))):
			# Convert the SAM file to a sorted BAM file.
			command = 'samtools view -bS %s | samtools sort - %s' % (alignments_path, os.path.splitext(alignments_path_bam)[0]);
			sys.stderr.write(command + '\n')
			subprocess.call(command, shell='True');
			# Create the BAM index file.
			command = 'samtools index %s %s.bai' % (alignments_path_bam, alignments_path_bam);
			subprocess.call(command, shell='True');		
	elif (alignments_path.endswith('bam') == False):
		sys.stderr.write('ERROR: File extension needs to be either .sam or .bam! Input file path: "%s".\n' % alignments_path);
		return;
	
	# Convert the sorted BAM file to a mpileup file if it doesn't exist yet.
	mpileup_path = ('%s.mpileup' % alignments_path_bam);
	mpileup_exists = os.path.exists(mpileup_path);
	if (mpileup_exists == False or (mpileup_exists == True and os.path.getmtime(alignments_path) > os.path.getmtime(mpileup_path))):
		command = 'samtools mpileup -B -d 1000000 -Q 0 -A -f %s %s > %s.mpileup' % (reference_path, alignments_path_bam, alignments_path_bam);
		subprocess.call(command, shell='True');

	sys.stderr.write('Processing file "%s"...\n' % alignments_path);
	sys.stderr.write('Coverage threshold: %d\n' % coverage_threshold);
	process_mpileup(alignments_path, reference_path, ('%s.mpileup' % alignments_path_bam), coverage_threshold, out_file, thread_id, bed_position);
	sys.stderr.write('\nMerging neighboring events...\n');
	chain_events(out_file);
	sys.stderr.write('Done!\n');

def CollectSummaries(sam_files, collective_output_file):
	fp_collect = None;
	
	try:
		fp_collect = open(collective_output_file, 'w');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % collective_output_file);
		return;
	
	for sam_file in sam_files:
		summary_file = os.path.splitext(sam_file)[0] + '.conssum';
		
		try:
			fp_sum = open(summary_file, 'r');
			lines = fp_sum.readlines();
			fp_sum.close();
		except IOError:
			sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % summary_file);
			continue;
		
		fp_collect.write(''.join(lines) + '\n');
	
	fp_collect.close();

if __name__ == "__main__":
	# if (len(sys.argv) < 5):
	# 	sys.stderr.write('Usage:\n');
	# 	sys.stderr.write('\t%s <reference_file_path> coverage_threshold <collective_output_file> <{sb}am_file_1> [<{sb}am_file_2> <{sb}am_file_3> ...]\n' % sys.argv[0]);
	# 	sys.stderr.write('\t(If <collective_output_file> is equal to "-", no files will be written to disk.)\n');
	# 	exit(1);

	if (len(sys.argv) < 4):
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <reference_file_path> <{sb}am_file> <out_file> [position]\n' % sys.argv[0]);
		sys.stderr.write('\tPosition parameter is a string specifying "chromosome:start-end"\n\n');
		exit(1);
	
	reference_file = sys.argv[1];
	coverage_threshold = 0; # int(sys.argv[2]);
	# output_prefix = sys.argv[3];
	output_prefix = '-';
	sam_file = sys.argv[2];
	out_file = sys.argv[3];
	bed_position = '';
	if (len(sys.argv) > 4):
		bed_position = sys.argv[4];
	# sys.stderr.write('bed_position: "%s"\n\n' % bed_position);
	
	processes = [];

	# if (output_prefix == '-'):
	# 	output_prefix = os.path.splitext(sam_file)[0];
	main(sam_file, reference_file, coverage_threshold, out_file, 0, bed_position);

	# if (output_prefix != '-'):
	# 	CollectSummaries([sam_file], output_prefix + '.variant.sum');
