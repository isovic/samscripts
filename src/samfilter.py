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

def filter_large_cigar_ops(sam_file, max_cigar_ops, out_filtered_sam_file):
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
		num_cigar_ops = sam_line.GetNumCigarOps();

		if (num_cigar_ops <= max_cigar_ops):
			fp_out.write(line);
			num_accepted += 1;

		else:
			sys.stderr.write('\tRejected: line: %d, qname: %s, len: %d, num_cigar_ops: %d, max_cigar_ops: %d\n' % (i, sam_line.qname, len(sam_line.seq), num_cigar_ops, max_cigar_ops));
			# print '\tRejected: line: %d, qname: %s, len: %d, CIGAR len: %d' % (i, sam_line.qname, len(sam_line.seq), cigar_len);
			num_rejected += 1;
		i += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('Number of good alignments: %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
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


def filter_qnames(sam_file, pattern_file, include_or_exclude, sam_parameter_to_compare, out_filtered_sam_file):
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

		if (include_or_exclude == 'include'):
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

		elif (include_or_exclude == 'exclude'):
			is_current_accepted = False;
			for pattern in patterns:
				if ((sam_parameter_to_compare == 'nqname' and pattern.lower() in sam_line.qname.lower()) or
					(sam_parameter_to_compare == 'nrname' and pattern.lower() in sam_line.rname.lower())):
					is_current_accepted = True;
					break;
			if (is_current_accepted == False):
				fp_out.write(line);
				num_accepted += 1;
				is_current_accepted = True;
			else:				
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

	try:
		sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
		sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));
	except:
		sys.stderr.write('num_accepted = %d (0.00 %%)\n' % (num_accepted));
		sys.stderr.write('num_rejected = %d (0.00 %%)\n' % (num_rejected));

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
		if ((is_1d_or_2d == '1d' and ((('twodir' in qname) == False and ('-2d' in qname) == False and ('_2d' in qname) == False) or ('-1d' in qname) == True or ('_1d' in qname) == True)) or
			(is_1d_or_2d == '2d' and (('twodir' in qname) == True or ('-2d' in qname) == True or ('_2d' in qname) == True))):

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


# def add_AS_from_NM(sam_file, out_filtered_sam_file):
# 	fp_in = None;
# 	fp_out = None;
	
# 	try:
# 		fp_in = open(sam_file, 'r');
# 	except IOError:
# 		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
# 		exit(1);
	
# 	try:
# 		fp_out = open(out_filtered_sam_file, 'w');
# 	except IOError:
# 		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_filtered_sam_file));
# 		exit(1);
	
# 	num_accepted = 0;
# 	num_rejected = 0;
	
# 	for line in fp_in:
# 		line = line.strip();
# 		if (len(line) == 0 or line[0] == '@'):
# 			fp_out.write(line);
# 			continue;
		
# 		if ('AS:i:' in line):
# 			fp_out.write(line);
# 			continue;
# 		if (('NM:i:' in line) == False):
# 			# fp_out.write(line);
# 			num_rejected += 1;
# 			continue;

# 		sam_line = utility_sam.SAMLine(line.rstrip());
# 		new_line = line + '\tAS:i:%d' % (sam_line.edit_distance);
# 		fp_out.write(new_line + '\n');

# 		# if (float(sam_line.mapq) >= float(mapq_limit)):
# 		# # if (float(sam_line.mapq) >= float(mapq_limit)):
# 		# 	fp_out.write(line);
# 		# 	num_accepted += 1;
# 		# else:
# 		# 	num_rejected += 1;
	
# 	fp_in.close();
# 	fp_out.close();
	
# 	# sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
# 	# sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def filter_out_of_bounds(sam_file, out_filtered_sam_file, out_rejected_sam_file=None):
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

	fp_out_rejected = None;
	try:
		fp_out_rejected = open(out_rejected_sam_file, 'w');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!' % (__name__, out_rejected_sam_file));
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
		
		sam_line = utility_sam.SAMLine(line.rstrip());

		pos_start = sam_line.pos - 1;
		len_on_ref = sam_line.CalcReferenceLengthFromCigar();
		pos_end = pos_start + len_on_ref - 1;
		reference_length = 4639675;

		# sam_line = utility_sam.SAMLine(line.rstrip());
		# split_line = line.rstrip().split('\t');
		# qname = split_line[0].lower();

		if (sam_line.IsMapped() == False or (sam_line.IsMapped() == True and pos_start >= 0 and pos_end < reference_length)):
			fp_out.write(line);
			num_accepted += 1;
		else:
			if (fp_out_rejected != None):
				fp_out_rejected.write(line);
			num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	if (fp_out_rejected != None):
		fp_out_rejected.close();
	
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def sam_stats(sam_file, reads_fastq=''):
	[hashed_sam, num_sam_lines, num_unique_sam_lines] = utility_sam.HashSAMWithFilter(sam_file, {});

	if (reads_fastq != ''):
		# [headers, seqs, quals] = fastqparser.read_fastq(reads_fastq);
		# header_hash = {};
		# i = 0;
		# while (i < len(headers)):
		# 	header_hash[headers[i]] = i;
		# 	header_hash[headers[i].split()[0]] = i;
		# 	i += 1;
		read_qnames = get_fastq_headers(reads_fastq);

	sys.stdout.write('Number of SAM lines in file: %d\n' % num_sam_lines);
	# sys.stdout.write('Number of SAM lines with unique qname: %d\n' % num_unique_sam_lines);

	num_mapped = 0;
	num_mapped_unique = 0;

	num_qnames_in_sam_file = len(hashed_sam.keys());

	for sam_lines in hashed_sam.values():
		for sam_line in sam_lines:
			if (sam_line.IsMapped() == True):
				num_mapped += 1;

		if (sam_lines[0].IsMapped() == True):
			num_mapped_unique += 1;

	sys.stdout.write('Mapped lines: %d / %d (%.2f%%)\n' % (num_mapped, num_sam_lines, float(num_mapped) / float(num_sam_lines) * 100.0));
	sys.stdout.write('Unique mapped lines: %d / %d (%.2f%%)\n' % (num_mapped_unique, num_qnames_in_sam_file, float(num_mapped_unique) / float(num_qnames_in_sam_file) * 100.0));
	if (reads_fastq != ''):
		sys.stdout.write('Mapped reads: %d / %d ( %.2f%%)\n' % (num_mapped_unique, len(read_qnames), float(num_mapped_unique) / float(len(read_qnames)) * 100.0));

def generate_AS(sam_file, reference_file, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty, out_filtered_sam_file):
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

	[ref_headers, ref_seqs, ref_quals] = fastqparser.read_fastq(reference_file);
	# [headers, seqs, quals] = fastqparser.read_fastq(fastq_file);
	ref_header_hash = {};
	i = 0;
	for header in ref_headers:
		ref_header_hash[header] = i;
		ref_header_hash[header.split(' ')[0]] = i;
		i += 1;

	num_accepted = 0;
	num_rejected = 0;
	num_unmapped = 0;
	
	i = 0;
	for line in fp_in:
		line = line.strip();
		if (len(line) == 0 or line[0] == '@' or (('AS:i:' in line) and ('NM:i:' in line))):
		# if (len(line) == 0 or line[0] == '@'):
			fp_out.write(line + '\n');
			continue;

		i += 1;
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		sam_line = utility_sam.SAMLine(line.rstrip());

		if (sam_line.IsMapped() == False):
			fp_out.write(line + '\n');
			num_unmapped += 1;
			continue;

		rname = sam_line.rname;
		try:
			ref_id = ref_header_hash[rname];
		except:
			sys.stderr.write('Warning: Could not find reference "%s"! Skipping alignment.\n' % (rname));
			fp_out.write(line + '\n');
			num_rejected += 1;
			continue;

		AS = 0;
		NM = 0;

		clip_start = sam_line.clip_count_front if (sam_line.clip_op_front == 'H') else 0;
		cigar_pos_list = sam_line.CalcCigarStartingPositions(separate_matches_in_individual_bases=True);
		for cigar_pos_info in cigar_pos_list:
			[cigar_count, cigar_op, pos_on_reference, pos_on_read] = cigar_pos_info;
			if (cigar_op == 'M'):
				if (sam_line.seq[pos_on_read - clip_start] == ref_seqs[ref_id][pos_on_reference]):
					AS += match_score;
				else:
					AS += mismatch_penalty;
					NM += 1;
			elif (cigar_op in 'ID'):
				AS += (gap_open_penalty + (cigar_count - 1) * gap_extend_penalty);
				NM += cigar_count;

		out_line = line;
		if (('AS:i:' in line) == False):
			out_line += '\tAS:i:%d' % (AS);
		if (('NM:i:' in line) == False):
			out_line += '\tNM:i:%d' % (NM);
		# if (sam_line.evalue <= threshold):
		# fp_out.write('%s\tAS:i:%d\tNM:i:%d\n' % (line, AS, NM));
		fp_out.write(out_line + '\n');
		num_accepted += 1;
		# else:
		# 	num_rejected += 1;
	
	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('Done!\n');
	try:
		sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
		sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));
		sys.stderr.write('num_unmapped = %d (%.2f%%)\n' % (num_unmapped, (float(num_unmapped) / float(num_accepted + num_rejected)) * 100.0));
	except:
		pass;

def marginalign_filter(sam_file, out_filtered_sam_file):
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
		line = line.strip();
		if (len(line) == 0 or line[0] == '@'):
			if (line.startswith('@SQ')):
				split_line = line.split('\t');
				found_hname = False;
				for param in split_line:
					if (param.startswith('SN:')):
						hname = param.split(':')[-1];
						new_hname = re.sub('[^0-9a-zA-Z]', '_', hname);
						sys.stderr.write('Found hname: "%s", replacing with: "%s".\n' % (hname, new_hname));
						new_line = line.replace(hname, new_hname);
						fp_out.write(new_line + '\n');
						found_hname = True;
						break;
				if (found_hname == False):
					fp_out.write(line + '\n');
			else:
				fp_out.write(line + '\n');
			continue;

		i += 1;
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		sam_line = utility_sam.SAMLine(line.rstrip());
		# split_line = line.split('\t');

		if (sam_line.IsMapped() == False):
			fp_out.write(line + '\n');
			num_unmapped += 1;

		qname = sam_line.qname;
		rname = sam_line.rname;

		new_line = line + '';

		if (qname != '*'):
			new_qname = re.sub('[^0-9a-zA-Z]', '_', qname);
			new_line = new_line.replace(qname, new_qname);

		if (rname != '*'):
			new_rname = re.sub('[^0-9a-zA-Z]', '_', rname);
			new_line = new_line.replace(rname, new_rname);

		fp_out.write(new_line + '\n');
		num_accepted += 1;

	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('Done!\n');
	try:
		sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
		sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));
		sys.stderr.write('num_unmapped = %d (%.2f%%)\n' % (num_unmapped, (float(num_unmapped) / float(num_accepted + num_rejected)) * 100.0));
	except:
		pass;

def fix_sam_hnames(sam_file, reference_path, out_filtered_sam_file):
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
		line = line.strip();
		if (len(line) == 0 or line[0] == '@'):
			if (line.startswith('@SQ')):
				split_line = line.split('\t');
				found_hname = False;
				for param in split_line:
					if (param.startswith('SN:')):
						hname = param.split(':')[-1];
						new_hname = hname.split()[0];
						sys.stderr.write('Found hname: "%s", replacing with: "%s".\n' % (hname, new_hname));
						new_line = line.replace(hname, new_hname);
						fp_out.write(new_line + '\n');
						found_hname = True;
						break;
				if (found_hname == False):
					fp_out.write(line + '\n');
			else:
				fp_out.write(line + '\n');
			continue;

		i += 1;
		sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		sam_line = utility_sam.SAMLine(line.rstrip());
		# split_line = line.split('\t');

		if (sam_line.IsMapped() == False):
			fp_out.write(line + '\n');
			num_unmapped += 1;

		qname = sam_line.qname;
		rname = sam_line.rname;

		new_line = line + '';

		if (qname != '*'):
#			new_qname = re.sub('[^0-9a-zA-Z]', '_', qname);
			new_qname = qname.split()[0];
			new_line = new_line.replace(qname, new_qname);

		if (rname != '*'):
#			new_rname = re.sub('[^0-9a-zA-Z]', '_', rname);
			new_rname = rname.split()[0];
			new_line = new_line.replace(rname, new_rname);

		fp_out.write(new_line + '\n');
		num_accepted += 1;

	fp_in.close();
	fp_out.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('Done!\n');
	try:
		sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
		sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));
		sys.stderr.write('num_unmapped = %d (%.2f%%)\n' % (num_unmapped, (float(num_unmapped) / float(num_accepted + num_rejected)) * 100.0));
	except:
		pass;

def sam_info(sam_file, reads_file=None):
	fp_in = None;
	fp_out = sys.stdout;
	
	try:
		fp_in = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, sam_file));
		exit(1);

	num_unmapped = 0;
	num_mapped_1d = 0;
	num_unmapped_1d = 0;
	num_mapped_2d = 0;
	num_unmapped_2d = 0;

	i = 0;
	for line in fp_in:
		line = line.strip();
		if (len(line) == 0 or line[0] == '@'):
			continue;

		i += 1;
		sys.stderr.write('\rLine %d' % (i));

		sam_line = utility_sam.SAMLine(line.rstrip());
		# split_line = line.split('\t');

		if (sam_line.IsMapped() == False):
			num_unmapped += 1;

		qname = sam_line.qname;
		rname = sam_line.rname;

		### Since there is no well-defined standard for naming the 1d and 2d reads, give precedance to
		### checking whether it's a template or a complement, and then check if it might be a 2d read.
		### If all this fails, then it probably is not a nanopore read, but from some other technology,
		### which means that it is definitely 1d.
		if (('template' in qname) or ('complement' in qname)):
			if (sam_line.IsMapped() == True):
				num_mapped_1d += 1;
			else:
				num_unmapped_1d += 1;
		elif (('2d' in qname) or ('2D' in qname) or ('twodir' in qname)):
			if (sam_line.IsMapped() == True):
				num_mapped_2d += 1;
			else:
				num_unmapped_2d += 1;
		else:
			if (sam_line.IsMapped() == True):
				num_mapped_1d += 1;
			else:
				num_unmapped_1d += 1;

	fp_in.close();

	sys.stderr.write('\n');
	fp_out.write('num_mapped_1d = %d\n' % (num_mapped_1d));
	fp_out.write('num_unmapped_1d = %d\n' % (num_unmapped_1d));
	fp_out.write('num_mapped_2d = %d\n' % (num_mapped_2d));
	fp_out.write('num_unmapped_2d = %d\n' % (num_unmapped_2d));

	try:
		fp_out.write('Percent of mapped reads which were 1d: %.2f%%\n' % (float(num_mapped_1d) / float(num_mapped_1d + num_mapped_2d) * 100.0));
		fp_out.write('Percent of mapped reads which were 2d: %.2f%%\n' % (float(num_mapped_2d) / float(num_mapped_1d + num_mapped_2d) * 100.0));
	except:
		pass;

	# fp_out.write('Percent mapped 2d over 1d: %.2f' % (float(num_mapped_2d) / float(num_mapped_1d)));
	# fp_out.write('Percent mapped 1d over 1d: %.2f' % (float(num_mapped_2d) / float(num_mapped_1d)));
	
	sys.stderr.write('\n');
	sys.stderr.write('Done!\n');
	try:
		sys.stderr.write('num_unmapped = %d (%.2f%%)\n' % (num_unmapped, (float(num_unmapped) / float(num_accepted + num_rejected)) * 100.0));
		# sys.stderr.write('num_mapped_1d = %d (%.2f%%)\n' % (num_unmapped, (float(num_unmapped) / float(num_accepted + num_rejected)) * 100.0));
	except:
		pass;

	sys.stderr.write('Counting mapped reads...\n');
	[num_alignments, num_mapped_alignments, num_unique_reads, num_mapped_reads, num_mapped_bases] = utility_sam.CountMappedReads(sam_file);
	if (reads_file != None):
		[fastqinfo_string, fastqinfo_num_seqs, fastqinfo_total_seq_len, fastqinfo_average_seq_len] = fastqparser.count_seq_length(reads_file);

	sys.stderr.write('\n');
	try:
		summary_read_count = '';
		summary_read_count += 'num_alignments: %d\n' % num_alignments;
		summary_read_count += 'num_mapped_alignments: %d (%.2f%%)\n' % (num_mapped_alignments, (float(num_mapped_alignments) / float(num_alignments)) * 100.0);
		summary_read_count += 'num_unmapped_alignments: %d (%.2f%%)\n' % ((num_alignments - num_mapped_alignments), (float(num_alignments - num_mapped_alignments) / float(num_alignments)) * 100.0);
		summary_read_count += 'num_mapped_reads: %d\n' % num_mapped_reads;
		summary_read_count += 'num_uniquely_mapped_reads: %d\n' % num_unique_reads;
		summary_read_count += 'num_mapped_bases: %d\n' % num_mapped_bases;

		if (reads_file != None):
			summary_read_count += 'num_read_in_input_reads_file: %d\n' % (fastqinfo_num_seqs);
			summary_read_count += 'num_bases_in_input_reads_file: %d\n' % (fastqinfo_total_seq_len);
			summary_read_count += 'percent_mapped_reads: %.2f%%\n' % ((float(num_mapped_reads) / float(fastqinfo_num_seqs)) * 100.0);
			summary_read_count += 'percent_mapped_bases: %.2f%%\n' % ((float(num_mapped_bases) / float(fastqinfo_total_seq_len)) * 100.0);
		fp_out.write(summary_read_count + '\n');
	except Exception, e:
		sys.stderr.write(str(e) + '\n');



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
		sys.stderr.write('\tnqname\n');
		sys.stderr.write('\tnrname\n');
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
		sys.stderr.write('\tstats\n');
		sys.stderr.write('\toutofbounds\n');
		sys.stderr.write('\tgenerateAS\n');
		sys.stderr.write('\tmarginalign\n');
		sys.stderr.write('\tfixhnames\n');
		sys.stderr.write('\tlongcigars\n');
		sys.stderr.write('\tinfo\n');


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
		filter_qnames(sam_file, patterns_file, 'include', sys.argv[1], out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'nqname' or sys.argv[1] == 'nrname'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Takes in a SAM file and a file containing qname patterns (one pattern in one line). Outputs all sam lines *not* containing a pattern.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_patterns_file> <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		patterns_file = sys.argv[2];
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_qnames(sam_file, patterns_file, 'exclude', sys.argv[1], out_filtered_sam_file);
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

		mapped = int(sys.argv[2]);
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

	elif (sys.argv[1] == 'outofbounds'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Check if alignment start and end positions fall within the boundaries of the reference genome.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file> <out_rejected_sam_file>\n' % (sys.argv[0], sys.argv[1], sys.argv[2]));
			exit(0);

		sam_file = sys.argv[2];
		out_filtered_sam_file = sys.argv[3];
		out_rejected_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file or sam_file == out_rejected_sam_file or out_filtered_sam_file == out_rejected_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		filter_out_of_bounds(sam_file, out_filtered_sam_file, out_rejected_sam_file);
		exit(0);

	elif (sys.argv[1] == 'generateAS'):
		if (len(sys.argv) != 5):
			sys.stderr.write('Calculates a simple alignment score if not previously present, using a model [1, -1, -1, -1].\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <reference_file> <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		reference_file = sys.argv[2];
		sam_file = sys.argv[3];
		out_filtered_sam_file = sys.argv[4];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		generate_AS(sam_file, reference_file, 1, -1, -1, -1, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'marginalign'):
		if (len(sys.argv) != 4):
			sys.stderr.write('Changes the qnames and the rnames of alignments not to include special characters.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		out_filtered_sam_file = sys.argv[3];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		marginalign_filter(sam_file, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'fixhnames'):
		if (len(sys.argv) < 4 or len(sys.argv) > 5):
			sys.stderr.write('Changes the qnames and the rnames of alignments not to include special characters.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file> [<reference_file>]\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		out_filtered_sam_file = sys.argv[3];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);
		reference_path = '';
		if (len(sys.argv) == 5):
			reference_path = sys.argv[5];
		fix_sam_hnames(sam_file, reference_path, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'longcigars'):
		if (len(sys.argv) != 4):
			sys.stderr.write('BAM file has a length limit on the number of CIGAR operations that it can store. The longer, erroneous reads might exceed these limits. Concretely, BAM format specification states:\n');
			sys.stderr.write('\tflag nc\tFLAG<<16|n cigar op;19 n cigar op is the number of operations\tuint32 t\n');
			sys.stderr.write('This means that the BAM file can hold up to 2^16-1 CIGAR operations.\n\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> <out_filtered_sam_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		out_filtered_sam_file = sys.argv[3];
		if (sam_file == out_filtered_sam_file):
			sys.stderr.write('ERROR: Output and input files are the same!\n');
			exit(0);

		filter_large_cigar_ops(sam_file, (2**16)-1, out_filtered_sam_file);
		exit(0);

	elif (sys.argv[1] == 'info'):
		if (len(sys.argv) < 3 or len(sys.argv) > 4):
			sys.stderr.write('Changes the qnames and the rnames of alignments not to include special characters.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_sam_file> [<reads_fastq_path>]\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		sam_file = sys.argv[2];
		reads_path = None;
		if (len(sys.argv) == 4):
			reads_path = sys.argv[3];
		sam_info(sam_file, reads_path);
		exit(0);



	else:
		sys.stderr.write('ERROR: Unknown subcommand!\n');
		exit(0);
