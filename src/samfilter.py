#! /usr/bin/python

import re;
import os;
import sys;
import utility_sam;
import fastqparser;



def filter_wrong_cigars(sam_file, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!\n' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	i = 0;
	for line in fp_in:
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));
		
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			i += 1;
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());
		# if (sam_line.line_fields_ok == False):
		# 	i += 1;
		# 	continue;

		# seq_len = len(sam_line.seq);
		# cigar_len = sam_line.CalcAlignmentLengthFromCigar();
		# if (sam_line.IsMapped() == False or (sam_line.IsMapped() == True and cigar_len == seq_len)):
		if (sam_line.IsAlignmentSane() == True):
			fp_out.write(line);
			num_accepted += 1;
		else:
			sys.stderr.write('\tRejected: line: %d, qname: %s, len: %d\n' % (i, sam_line.qname, len(sam_line.seq)));
			# print '\tRejected: line: %d, qname: %s, len: %d, CIGAR len: %d' % (i, sam_line.qname, len(sam_line.seq), cigar_len);
			num_rejected += 1;
		i += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('Number of alignments with equal seq length and CIGAR length: %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('Number of faulty alignments: %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def filter_sam_by_mapq(sam_file, mapq_limit, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());
		
		if (float(sam_line.mapq) >= float(mapq_limit)):
		# if (float(sam_line.mapq) >= float(mapq_limit)):
			fp_out.write(line);
			num_accepted += 1;
		else:
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def filter_sam_by_as(sam_file, as_limit, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		alignment_score = '';
		try:
			alignment_score = sam_line.optional['AS'];
		except:
			continue;
		
		if (float(alignment_score) >= float(as_limit)):
			fp_out.write(line);
			num_accepted += 1;
		else:
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def filter_sam_by_length(sam_file, threshold, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		if ((len(sam_line.seq) + sam_line.clip_count_front + sam_line.clip_count_back) >= threshold):
			fp_out.write(line);
			num_accepted += 1;
		else:
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def filter_sam_by_evalue(sam_file, threshold, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	num_unmapped = 0;
	
	i = 0;
	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;

		i += 1;
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		# if ((len(sam_line.seq) + sam_line.clip_count_front + sam_line.clip_count_back) >= threshold):
		if (sam_line.IsMapped() == False):
			num_unmapped += 1;

		if (sam_line.evalue <= threshold):
			fp_out.write(line);
			num_accepted += 1;
		else:
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('Done!\n');
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_unmapped = %d (%.2f%%)\n' % (num_unmapped, (float(num_unmapped) / float(num_accepted + num_rejected)) * 100.0));


def filter_qnames(sam_file, pattern_file, sam_parameter_to_compare, out_filtered_sam_file):
	fp_in = None;
	fp_patterns = None;
	fp_out = None;

	try:
		fp_patterns = open(pattern_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, pattern_file));
		exit(1);

	patterns = [line.strip() for line in fp_patterns.readlines() if (len(line.strip()) > 0)];
	fp_patterns.close();

	if (len(patterns) == 0):
		sys.stderr.write('No qname patterns given. Exiting.\n');
		return;

	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		is_current_accepted = False;
		for pattern in patterns:
			if ((sam_parameter_to_compare == 'qname' and pattern.lower() in sam_line.qname.lower()) or
				(sam_parameter_to_compare == 'rname' and pattern.lower() in sam_line.rname.lower())):

				fp_out.write(line);
				num_accepted += 1;
				is_current_accepted = True;
				break;
			# else:
				# print 'Tu sam 1!';
				# print pattern;
				# print sam_line.qname;
				# print sam_line.rname;
		if (is_current_accepted == False):
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

# def filter_uniqe_best(sam_file, out_filtered_sam_file):
# 	[headers, sam_lines] = utility_sam.LoadSAM(sam_file);

# 	fp_out = None;
# 	try:
# 		fp_out = open(out_filtered_sam_file, 'w');
# 	except IOError:
# 		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_filtered_sam_file));
# 		exit(1);

	
# 	sorted_lines_by_name = sorted(sam_lines, key=lambda sam_line: (sam_line.qname, (-sam_line.chosen_quality)));

# 	# Filter unique SAM lines, where uniqueness is defined by the qname parameter.
# 	# If there is more than one alignment with the same name, pick only the first one
# 	# because we have already sorted them by quality.
# 	unique_lines = [];
# 	i = 0;
# 	previous_line = None;
# 	num_unambiguous_reads = 0;
# 	num_unmapped_alignments = 0;
# 	num_unmapped_reads = 0;
# 	num_unambiguous_mapped_reads = 0;
# 	is_read_mapped = 0;
# 	count = 0;
# 	while i < len(sorted_lines_by_name):
# 		sam_line = sorted_lines_by_name[i];
		
# 		if (sam_line.IsMapped() == False):
# 			num_unmapped_alignments += 1;
# 		else:
# 			is_read_mapped = 1;
		
# 		if (previous_line == None or sam_line.qname != previous_line.qname):
# 			if (count == 1):
# 				num_unambiguous_reads += 1;
# 				num_unambiguous_mapped_reads += is_read_mapped;
# 			unique_lines.append(sam_line);
# 			num_unmapped_reads += 1 if (sam_line.IsMapped() == False) else 0;
# 			count = 1;
# 			is_read_mapped = 0;
# 		elif (sam_line.qname == previous_line.qname):
# 			count += 1;
# 		previous_line = sam_line;
# 		i += 1;
# 	if (count == 1):
# 		num_unambiguous_reads += 1;
# 		num_unambiguous_mapped_reads += is_read_mapped;

# 	fp_out.write('\n'.join(headers) + '\n');
# 	fp_out.write('\n'.join([sam_line.original_line for sam_line in unique_lines]));

# 	num_accepted = len(unique_lines);
# 	num_rejected = len(sam_lines) - len(unique_lines);

# 	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
# 	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

# 	fp_out.close();

def filter_uniqe_best(sam_file, out_filtered_sam_file, out_rejected_file=None):

	fp_out = None;
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_filtered_sam_file));
		exit(1);

	fp_out_rejected = None;
	try:
		if (out_rejected_file != None):
			fp_out_rejected = open(out_rejected_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_rejected_file));
		exit(1);

	[hashed_sam, num_sam_lines, num_unique_sam_lines] = utility_sam.HashSAMWithFilter(sam_file, {});
	sam_headers = utility_sam.LoadOnlySAMHeaders(sam_file);

	fp_out.write('\n'.join(sam_headers) + '\n');

	num_accepted = 0;
	num_rejected = 0;

	i = 0;
	for sam_lines in hashed_sam.values():
		i += 1;
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));
		sam_line = sam_lines[0];
		fp_out.write(sam_line.original_line + '\n');
		num_accepted += 1;
		num_rejected += len(sam_lines) - 1;

		if (fp_out_rejected != None):
			for sam_line in sam_lines[1:]:
				fp_out_rejected.write(sam_line.original_line + '\n');



	sys.stderr.write('\n');

	# num_accepted = len(hashed_sam.keys());
	# num_rejected = num_sam_lines - num_accepted;

	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

	if (fp_out != None):
		fp_out.close();
	if (fp_out_rejected != None):
		fp_out_rejected.close();

def filter_uniqe_best_with_AS_0(sam_file, out_filtered_sam_file):

	fp_out = None;
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_filtered_sam_file));
		exit(1);

	[hashed_sam, num_sam_lines, num_unique_sam_lines] = utility_sam.HashSAMWithFilter(sam_file, {});
	sam_headers = utility_sam.LoadOnlySAMHeaders(sam_file);

	fp_out.write('\n'.join(sam_headers) + '\n');

	for sam_lines in hashed_sam.values():
		if (len(sam_lines) == 1):
			sam_line = sam_lines[0];
			fp_out.write(sam_line.original_line + '\n');
		else:
			if (sam_lines[0].AS > 0):
				sam_line = sam_lines[0];
				fp_out.write(sam_line.original_line + '\n');

	num_accepted = len(hashed_sam.keys());
	num_rejected = num_sam_lines - num_accepted;

	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

	fp_out.close();



def filter_mapped(sam_file, mapped, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		if ((mapped == 0 and sam_line.IsMapped() == False) or (mapped != 0 and sam_line.IsMapped() == True)):
			fp_out.write(line);
			num_accepted += 1;
		else:
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def filter_mapped_qnames(sam_file, mapped, out_filtered_sam_file=None):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);

	fp_out = sys.stdout;
	if (out_filtered_sam_file != None):
		try:
			fp_out = open(out_filtered_sam_file, 'w');
		except IOError:
			sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
			exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			# fp_out.write(line);
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		if ((mapped == 0 and sam_line.IsMapped() == False) or (mapped != 0 and sam_line.IsMapped() == True)):
			# fp_out.write(line);
			fp_out.write(sam_line.qname + '\n');
			num_accepted += 1;
		else:
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def add_quality_values(sam_file, fastq_file, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);

	[headers, seqs, quals] = fastqparser.read_fastq(fastq_file);
	header_hash = {};
	i = 0;
	while (i < len(headers)):
		header_hash[headers[i]] = i;
		header_hash[headers[i].split()[0]] = i;
		i += 1;

	num_accepted = 0;
	num_rejected = 0;
	
	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		if (sam_line.qual != '*' and sam_line.qual != ''):
			fp_out.write(line);
			continue;

		read_id = header_hash[sam_line.qname];
		read_qual = quals[read_id];
		if (read_qual == None or len(read_qual) == 0):
			sys.stderr.write('ERROR: Reads file given is not FASTQ! Exiting.\n');
			exit(1);

		clip_start = 0;
		if (sam_line.clip_op_front == 'H'):
			clip_start = sam_line.clip_count_front;
		clip_end = len(read_qual);
		if (sam_line.clip_op_back == 'H'):
			clip_end -= sam_line.clip_count_back;

		read_clipped_qual = read_qual[clip_start:clip_end];
		split_original_line = sam_line.original_line.strip().split('\t');
		split_original_line[10] = read_clipped_qual;
		new_line = '\t'.join(split_original_line);
		fp_out.write(new_line + '\n');

		# print sam_line.cigar;
		# print read_qual;
		# print read_qual[0:(clip_start)], ' ', read_clipped_qual, '\n', read_qual[clip_end:];
		# exit(1);

	fp_in.close();
	fp_out.close();

def add_dummy_quality_values(phred_qv, sam_file, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);

	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		split_line = line.strip().split('\t');
		seq = split_line[9];
		split_line[10] = phred_qv * len(seq);
		new_line = '\t'.join(split_line);
		fp_out.write(new_line + '\n');

	fp_in.close();
	fp_out.close();

def add_dummy_quality_values_nanopore(phred_qv_1d, phred_qv_2d, sam_file, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		split_line = line.strip().split('\t');
		qname = split_line[0].lower();
		seq = split_line[9];
		if (('twodir' in qname) or ('-2d' in qname)):
			split_line[10] = phred_qv_2d * len(seq);
		else:
			split_line[10] = phred_qv_1d * len(seq);
		new_line = '\t'.join(split_line);
		fp_out.write(new_line + '\n');

	fp_in.close();
	fp_out.close();



def extract_region_full(region, sam_file, out_filtered_sam_file):
	try:
		[chromosome, location] = region.split(':');
		[start, end] = location.split('-');
		start = int(start);
		end = int(end);
	except:
		sys.stderr.write('ERROR: Region not formatted correctly! Given region: "%s", correct format: "chromosome:start-end.\n', region);
		exit(1);

	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	i = 0;
	for line in fp_in:
		i += 1;
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		# if ((mapped == 0 and sam_line.IsMapped() == False) or (mapped != 0 and sam_line.IsMapped() == True)):
		if (sam_line.IsMapped() == True and sam_line.rname.startswith(chromosome) == True and sam_line.pos <= start and (sam_line.pos + sam_line.CalcReferenceLengthFromCigar()) >= end):
			fp_out.write(line);
			num_accepted += 1;
		else:
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def extract_region_partial(region, sam_file, out_filtered_sam_file):
	try:
		[chromosome, location] = region.split(':');
		[start, end] = location.split('-');
		start = int(start);
		end = int(end);
	except:
		sys.stderr.write('ERROR: Region not formatted correctly! Given region: "%s", correct format: "chromosome:start-end.\n', region);
		exit(1);

	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	i = 0;
	for line in fp_in:
		i += 1;
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		# if ((mapped == 0 and sam_line.IsMapped() == False) or (mapped != 0 and sam_line.IsMapped() == True)):
		alignment_start = sam_line.pos;
		alignment_end = sam_line.pos + sam_line.CalcReferenceLengthFromCigar();
		if (sam_line.IsMapped() == False or sam_line.rname.startswith(chromosome) == False or
			(alignment_start < start and alignment_end < start) or (alignment_start > end and alignment_end > end)):
			num_rejected += 1;
		else:
			fp_out.write(line);
			num_accepted += 1;

		# if (sam_line.IsMapped() == True and sam_line.rname.startswith(chromosome) == True and
		# 	( (alignment_start <= start and alignment_end >= start) or
		# 	  (alignment_start )):

		# 	fp_out.write(line);
		# 	num_accepted += 1;
		# else:
		# 	num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));
	


def filter_1d_2d(is_1d_or_2d, sam_file, out_filtered_sam_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_filtered_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;
	
	i = 0;
	for line in fp_in:
		i += 1;
		if (len(line.strip()) == 0 or line[0] == '@'):
			fp_out.write(line);
			continue;
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));
		
		# sam_line = utility_sam.SAMLine(line.rstrip());
		split_line = line.rstrip().split('\t');
		qname = split_line[0].lower();
		if ((is_1d_or_2d == '1d' and ((('twodir' in qname) == False and ('-2d' in qname) == False) or ('-1d' in qname) == True)) or
			(is_1d_or_2d == '2d' and (('twodir' in qname) == True or ('-2d' in qname) == True))):

			fp_out.write(line);
			num_accepted += 1;
		else:
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));



# /home/isovic/work/eclipse-workspace/data/graphmap-params_test-drive.sam /home/isovic/work/eclipse-workspace/data/graphmap-params_test-drive.csv
def extract_features_csv(sam_file, out_csv_file):
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);
	
	try:
		fp_out = open(out_csv_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_csv_file));
		exit(1);

	feature_labels = [];
	all_feature_lines = [];

	for line in fp_in:
		if (len(line.strip()) == 0 or line[0] == '@'):
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		if (sam_line.IsMapped() == False):
			continue;

		feature_line = [];
		feature_labels = [];

		feature_line.append((sam_line.qname));
		feature_labels.append('qname');

		feature_line.append('1' if ('escherichia coli' in sam_line.rname.lower()) else '0')
		feature_labels.append('is_good');

		feature_line.append((sam_line.rname));
		feature_labels.append('rname');

		# feature_line.append(str(sam_line.mapq));
		# feature_labels.append('mapq');

		feature_line.append(str(len(sam_line.seq)));
		feature_labels.append('seq_len');

		# feature_line.append(str(sam_line.pos));
		# feature_labels.append('position');

		# if ('gi|161501984|ref|NC_010067.1| Salmonella enterica subsp. arizonae serovar 62:z4' in sam_line.rname):
		# 	print sam_line.Verbose();
			# print sam_line.original_line;

		# split_line = sam_line.original_line.strip().split();

		# parameter = 'AS_final';
		# if (parameter in sam_line.optional):
		# 	feature_line.append(str(sam_line.optional[parameter]));
		# 	feature_labels.append(parameter);

		parameter = 'H0';
		if (parameter in sam_line.optional):
			feature_line.append(str(sam_line.optional[parameter]));
			feature_labels.append(parameter);

		# parameter = 'X9';
		# if (parameter in sam_line.optional):
		# 	feature_line.append(str(sam_line.optional[parameter]));
		# 	feature_labels.append(parameter);

		parameter = 'NM';
		NM = 0;
		if (parameter in sam_line.optional):
			feature_line.append(str(sam_line.optional[parameter]));
			feature_labels.append(parameter);
			NM = sam_line.optional[parameter];

		parameter = 'Y1';
		read_length = 0;
		if (parameter in sam_line.optional):
			read_length = sam_line.optional[parameter].split('_')[-1]
			feature_line.append(str(read_length));
			feature_labels.append('read_length');

		if (read_length != 0):
			feature_line.append(str(float(NM) / (float(read_length) / 2.0)));
			feature_labels.append('NM_ratio');

		parameter = 'X1';
		if (parameter in sam_line.optional):
			X1 = sam_line.optional[parameter];

			subparam = 'lcs_length';
			m1 = re.match('.*' + subparam + '=([0-9]+).*', X1);
			if (m1 != None):
				feature_line.append(str(m1.groups()[0]));
				feature_labels.append(subparam);

			subparam = 'cov_bases';
			m1 = re.match('.*' + subparam + '=([0-9]+).*', X1);
			if (m1 != None):
				feature_line.append(str(m1.groups()[0]));
				feature_labels.append(subparam);

			subparam = 'act_q';
			m1 = re.match('.*' + subparam + '\[([0-9]+,[0-9]+)\].*', X1);
			if (m1 != None):
				query = [int(value) for value in m1.groups()[0].split(',')];
				# feature_line.append(str((query[1] - query[0])));
				# feature_labels.append(act_query_length);

			subparam = 'act_r';
			m1 = re.match('.*' + subparam + '\[([0-9]+,[0-9]+)\].*', X1);
			if (m1 != None):
				reference = [int(value) for value in m1.groups()[0].split(',')];
				# feature_line.append(str((query[1] - query[0])));
				# feature_labels.append(act_query_length);

			subparam = 'act_d';
			m1 = re.match('.*' + subparam + '\[([0-9]+,[0-9]+)\].*', X1);
			if (m1 != None):
				distance = [int(value) for value in m1.groups()[0].split(',')];
				feature_line.append(str(distance[0]));
				feature_labels.append('d_query');
				feature_line.append(str(distance[1]));
				feature_labels.append('d_reference');

			# subparam = 'p';
			# m1 = re.match('.*_' + subparam + '\[_([0-9\.]+)\].*', X1);
			# if (m1 != None):
			# 	feature_line.append(str(m1.groups()[0]));
			# 	feature_labels.append(subparam);

			subparam = 'supp';
			m1 = re.match('.*' + subparam + '\[([0-9\.]+)\].*', X1);
			if (m1 != None):
				feature_line.append(str(m1.groups()[0]));
				feature_labels.append(subparam);

			subparam = 'std';
			m1 = re.match('.*\]_' + subparam + '\[([0-9\.]+)\].*', X1);
			if (m1 != None):
				feature_line.append(str(m1.groups()[0]));
				feature_labels.append(subparam);

			# subparam = 'AS';
			# m1 = re.match('.*' + subparam + '\[([0-9\.]+)\].*', X1);
			# if (m1 != None):
			# 	feature_line.append(str(m1.groups()[0]));
			# 	feature_labels.append(subparam);

			# subparam = 'AS_std';
			# m1 = re.match('.*' + subparam + '\[([0-9\.]+)\].*', X1);
			# if (m1 != None):
			# 	feature_line.append(str(m1.groups()[0]));
			# 	feature_labels.append(subparam);

			# subparam = 'AS_cov_bases';
			# m1 = re.match('.*' + subparam + '\[([0-9\.]+)\].*', X1);
			# if (m1 != None):
			# 	feature_line.append(str(m1.groups()[0]));
			# 	feature_labels.append(subparam);

			# subparam = 'AS_read_len';
			# m1 = re.match('.*' + subparam + '\[([0-9\.]+)\].*', X1);
			# if (m1 != None):
			# 	feature_line.append(str(m1.groups()[0]));
			# 	feature_labels.append(subparam);

			# subparam = 'AS_query_len';
			# m1 = re.match('.*' + subparam + '\[([0-9\.]+)\].*', X1);
			# if (m1 != None):
			# 	feature_line.append(str(m1.groups()[0]));
			# 	feature_labels.append(subparam);

			# subparam = 'edit_distance';
			# m1 = re.match('.*__' + subparam + '_([0-9]+).*', X1);
			# if (m1 != None):
			# 	feature_line.append(str(m1.groups()[0]));
			# 	feature_labels.append(subparam);

		all_feature_lines.append('\t'.join(feature_line));

	fp_out.write('\t'.join(feature_labels) + '\n');
	fp_out.write('\n'.join(all_feature_lines));

	fp_in.close();
	fp_out.close();

# AS:i:56 H0:i:2  X9:i:1  NM:i:4481       Y1:Z:read_length_12373  X1:Z:_best_path__lcs_length=5659_cov_bases=4051_local_scores_id_5_act_q[10,12299]_act_r[787651203,787664587]_act_d[12289,13384]_p[_0.918186]
# _supp[0.0818141]_std[220.923]_AS[0.5619]_AS_std[0.9439]_AS_cov_bases[0.5994]_AS_read_len[1.0000]_AS_query_len[0.9932]_mean_edit_distance_5567.85__max_edit_distance_6805.15__edit_distance_4481_


if __name__ == "__main__":
	if (len(sys.argv) < 2):
		sys.stderr.write('Various filtering methods for SAM files.\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\tmapq\n');
		sys.stderr.write('\tAS\n');
		sys.stderr.write('\tlength\n');
		sys.stderr.write('\twrongcigars\n');
		sys.stderr.write('\tqname\n');
		sys.stderr.write('\trname\n');
		sys.stderr.write('\tuniquebest\n');
		sys.stderr.write('\tmapped\n');
		sys.stderr.write('\tmappedqnames\n');
		sys.stderr.write('\tfeaturescsv\n');
		sys.stderr.write('\taddqv\n');
		sys.stderr.write('\tadddummyqv\n');
		sys.stderr.write('\tadddummyqvnanopore\n');
		sys.stderr.write('\tregionfull\n');
		sys.stderr.write('\tregionpartial\n');
		sys.stderr.write('\tevalue\n');
		sys.stderr.write('\t1d\n');
		sys.stderr.write('\t2d\n');
		exit(0);

	if (sys.argv[1] == 'mapq'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Extract SAM lines which have mapping quality higher or equal to a given threshold.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s mapq_limit <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		mapq_limit = float(sys.argv[2]);
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_sam_by_mapq(sam_file, mapq_limit, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'AS'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Extract SAM lines which have alignment score (AS field) higher or equal to a given threshold.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s as_limit <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		as_limit = float(sys.argv[2]);
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_sam_by_as(sam_file, as_limit, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'length'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Extract SAM lines which have sequence length higher or equal to a given threshold.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s len_limit <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		length_limit = float(sys.argv[2]);
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_sam_by_length(sam_file, length_limit, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'wrongcigars'):
		if (len(sys.argv) != 4):
			sys.stderr.write('Takes an input SAM file and for each entry compares the SEQ length and the length obtained from counting the CIGAR operations. These should match. If they do not match, these lines are filtered out. The output SAM contains only lines with matching lengths.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		out_filtered_sam_file = sys.argv[3];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_wrong_cigars(sam_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'qname' or sys.argv[1] == 'rname'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Takes in a SAM file and a file containing qname patterns (one pattern in one line). Outputs all sam lines containing a pattern.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_patterns_file> <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		patterns_file = sys.argv[2];
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_qnames(sam_file, patterns_file, sys.argv[1], out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'uniquebest'):
		if (len(sys.argv) < 4 or len(sys.argv) > 5):
			sys.stderr.write('Filteres only unique qname lines. If there are multiple alignments with same qname, then only the one with highest mapq/alignment score will be output.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		out_filtered_sam_file = sys.argv[3];
		out_rejected_file = None;
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		if (len(sys.argv) == 5):
			out_rejected_file = sys.argv[4];
			if (sam_file == out_rejected_file):
				sys.stderr.write('ERROR: Output for rejected SAM lines and input files are the same!\n');
				exit(0);
		filter_uniqe_best(sam_file, out_filtered_sam_file, out_rejected_file);
		exit(0);

	elif (sys.argv[1] == 'uniquebestas0'):
		if (len(sys.argv) != 4):
			sys.stderr.write('Filteres only unique qname lines. If there are multiple alignments with same qname, the ones with AS <= 0 are discarded, and then only the one with highest mapq/alignment score will be output.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		out_filtered_sam_file = sys.argv[3];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_uniqe_best_with_AS_0(sam_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'mapped'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Output only mapped (or unmapped) alignments. If parameter is 0, unmapped will be output, otherwise only mapped.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s mapped <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		mapped = sys.argv[2];
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_mapped(sam_file, mapped, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'mappedqnames'):
		if (len(sys.argv) < 4 or len(sys.argv) > 5):
			sys.stderr.write('Output only qnames of mapped (or unmapped) alignments. If parameter is 0, unmapped will be output, otherwise only mapped.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s mapped <input_sam_file> [<out_filtered_sam_file>]\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		mapped = sys.argv[2];
		sam_file = sys.argv[3];
		out_filtered_sam_file = None;
		if (len(sys.argv) == 5):
			out_filtered_sam_file = sys.argv[4];
			if (sam_file == out_filtered_sam_file):
				sys.stderr.write('ERROR: Output and input files are the same!\n');
				exit(0);
		filter_mapped_qnames(sam_file, mapped, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == '1d'):
		if (len(sys.argv) != 4):
			sys.stderr.write('Extracts only 1d reads from the given SAM file. If the qname does not contain keywords "twodir" or "-2d" or contains "-1d", it is considered 1d.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		out_filtered_sam_file = sys.argv[3];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_1d_2d(sys.argv[1], sam_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == '2d'):
		if (len(sys.argv) != 4):
			sys.stderr.write('Extracts only 2d reads from the given SAM file. If the qname contains keywords "twodir" or "-2d" it is considered as 2d.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		out_filtered_sam_file = sys.argv[3];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_1d_2d(sys.argv[1], sam_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'featurescsv'):
		if (len(sys.argv) != 4):
			sys.stderr.write('Extract various features from the SAM file alignments, and create a CSV formatted output suitable for analysis with RapidMiner.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		out_csv_file = sys.argv[3];
		if (sam_file == out_csv_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		extract_features_csv(sam_file, out_csv_file);
		exit(0);

	elif (sys.argv[1] == 'addqv'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Adds base quality values to the SAM lines from the given reads FASTQ file.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <reads_fastq> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		fastq_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		add_quality_values(sam_file, fastq_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'adddummyqv'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Adds base quality values to the SAM lines from the given reads FASTQ file.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s phredqv <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			sys.stderr.write('\n');
			sys.stderr.write('\t\tphredqv - decimal value of the Phred substitution score\n');
			exit(0);

		phred_qv = chr(int(sys.argv[2]) + 33);
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		add_dummy_quality_values(phred_qv, sam_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'adddummyqvnanopore'):
		if (len(sys.argv) != 6):
			sys.stderr.write('Adds base quality values to the SAM lines from the given reads FASTQ file.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s phredqv_default phredqv_1d phredqv_2d <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			sys.stderr.write('\n');
			sys.stderr.write('\t\tphredqv_1d - average decimal value of the Phred substitution score for 1d reads [18%% is 7.447 Phred]\n');
			sys.stderr.write('\t\tphredqv_2d - average decimal value of the Phred substitution score for 2d reads [11%% is 9.586 Phred]\n');
			exit(0);

		phred_qv_1d = chr(int(sys.argv[2]) + 33);
		phred_qv_2d = chr(int(sys.argv[3]) + 33);
		sam_file = sys.argv[4];
		out_filtered_sam_file = sys.argv[5];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		add_dummy_quality_values_nanopore(phred_qv_1d, phred_qv_2d, sam_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'regionfull'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Extracts only alignments that are fully overlapping given genomic region (in SAMtools format, i.e. chrom_name:start-end).\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s region <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		region =sys.argv[2];
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		extract_region_full(region, sam_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'regionpartial'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Extracts only alignments that are fully overlapping given genomic region (in SAMtools format, i.e. chrom_name:start-end).\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s region <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		region =sys.argv[2];
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		extract_region_partial(region, sam_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'evalue'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Extract SAM lines which have E-value less than a given threshold (custom field ZE in SAM output).\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s threshold <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		evalue_threshold = float(sys.argv[2]);
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_sam_by_evalue(sam_file, evalue_threshold, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'stats'):
		if (len(sys.argv) != 3):
			sys.stderr.write('Calculates some simple statistics on a SAM file, such as number of unique alignment, number of mapped alignments, and similar.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		sam_stats(sam_file);
		exit(0);

	else:
		sys.stderr.write('ERROR: Unknown subcommand!\n');
		exit(0);
