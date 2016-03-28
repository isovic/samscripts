#! /usr/bin/python

# Written by Ivan Sovic, March 2015.
#
# The MIT License (MIT)

# Copyright (c) 2015 Ivan Sovic

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



import os;
import sys;
import utility_sam;

def CompareTwoSAMs(sam_file1, sam_file2, distance_threshold, out_summary_prefix=''):

	# print 'Loading first SAM file...';
	qnames_with_multiple_alignments = {};
	# print sam_file1, sam_file2, distance_threshold, out_summary_prefix
	sys.stderr.write('Loading the first SAM file into hash...\n');
	[sam_hash1, sam_hash1_num_lines, sam_hash1_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_file1, qnames_with_multiple_alignments);
	sam_headers1 = utility_sam.LoadOnlySAMHeaders(sam_file1);
	sys.stderr.write('Loading the second SAM file into hash...\n');
	[sam_hash2, sam_hash2_num_lines, sam_hash2_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_file2, qnames_with_multiple_alignments);
	sam_headers2 = utility_sam.LoadOnlySAMHeaders(sam_file2);

	not_in_sam_file1 = 0;
	not_in_sam_file2 = 0;

	num_different_reference = 0;
	num_different_orientation = 0;
	num_not_mapped_1 = 0;
	num_not_mapped_2 = 0;
	num_mapped_1 = 0;
	num_mapped_2 = 0;

	qname_to_distance_hash = {};
	# qname_to_pos = {};
	distance_count_hash = {};
	distance_to_qname_hash = {};
	distance_to_sam_hash = {};
	shared_qnames = {};

	num_processed = 0;

	qnames_not_in_sam_file1 = [];
	qnames_not_in_sam_file2 = [];
	qnames_different_ref = [];
	qnames_different_orient = [];

	for qname in sam_hash1.keys():
		num_processed += 1;
		if ((num_processed % 1000) == 0):
			sys.stderr.write('\rProcessed %d alignments...' % num_processed);

		if (len(sam_hash1[qname]) > 0 and sam_hash1[qname][0].IsMapped() == True):
			num_mapped_1 += 1;

		# TODO: THIS NEEDS TO BE REMOVED OR IMPLEMENTED SOMEHOW DIFFERENTLY!!
		# The point of this was that, BLASR doesn't conform to the SAM standard, and makes it difficult to
		# uniformly evaluate the results!
		# if 'blasr' in sam_file1.lower():
		# 	qname = '/'.join(qname.split('/')[:-1]);

		sam_line_list1 = sam_hash1[qname];

		sam_line_list_2 = [];
		try:
			sam_line_list2 = sam_hash2[qname];
			if (len(sam_line_list2) == 0 or (len(sam_line_list2) > 0 and sam_line_list2[0].IsMapped() == False)):
				not_in_sam_file2 += 1;
				qnames_not_in_sam_file2.append([sam_line_list1[0].evalue, qname]);
				continue;
		except:
			not_in_sam_file2 += 1;
			qnames_not_in_sam_file2.append([sam_line_list1[0].evalue, qname]);
			continue;

		sorted_sam_line_list1 = sorted(sam_line_list1, key=lambda sam_line: ((not sam_line.IsMapped()), -sam_line.chosen_quality));
		sorted_sam_line_list2 = sorted(sam_line_list2, key=lambda sam_line: ((not sam_line.IsMapped()), -sam_line.chosen_quality));
		if (len(sorted_sam_line_list1) > 0 and len(sorted_sam_line_list2) > 0):

			if (sorted_sam_line_list1[0].IsMapped() == False):
				num_not_mapped_1 += 1;
			if (sorted_sam_line_list2[0].IsMapped() == False):
				num_not_mapped_2 += 1;
			if (sorted_sam_line_list1[0].IsMapped() == False or sorted_sam_line_list2[0].IsMapped() == False):
				continue;

			if (not ((sorted_sam_line_list1[0].rname in sorted_sam_line_list2[0].rname) or (sorted_sam_line_list2[0].rname in sorted_sam_line_list1[0].rname))):
				num_different_reference += 1;
				qnames_different_ref.append(qname);
				continue;
			if (sorted_sam_line_list1[0].IsReverse() != sorted_sam_line_list2[0].IsReverse()):
				num_different_orientation += 1;
				qnames_different_orient.append(qname);
				continue;

			distance = abs(sorted_sam_line_list1[0].clipped_pos - sorted_sam_line_list2[0].clipped_pos);
			if (not (qname in shared_qnames)):
				shared_qnames[qname] = 1;
			qname_to_distance_hash[qname] = distance;
			# qname_to_pos[qname] = [sorted_sam_line_list1[0].clipped_pos, sorted_sam_line_list2[0].clipped_pos];
			if (distance in distance_count_hash):
				distance_count_hash[distance] += 1;
				distance_to_qname_hash[distance].append(qname);
				distance_to_sam_hash[distance].append(sorted_sam_line_list1[0]);
			else:
				distance_count_hash[distance] = 1;
				distance_to_qname_hash[distance] = [qname];
				distance_to_sam_hash[distance] = [sorted_sam_line_list1[0]];
		else:
			if (len(sorted_sam_line_list1) == 0):
				not_in_sam_file1 += 1;
				# qnames_not_in_sam_file1.append(qname);
			if (len(sorted_sam_line_list2) == 0):
				not_in_sam_file2 += 1;
				# qnames_not_in_sam_file2.append(qname);
			sys.stderr.write('Warning: Something odd with qname "%s". Qname present in both files, but lists are empty.\n' % (qname));
			continue;

		# min_distance = -1;
		# i = 0;
		# for sam_line1 in sorted_sam_line_list1:
		# 	for sam_line2 in sorted_sam_line_list2:
		# 		distance = abs(sam_line1.clipped_pos - sam_line2.clipped_pos);
		# 		if (i == 0 or distance < min_distance):
		# 			min_distance = distance;
		# 		i += 1;
		# distance_hash[qname] = min_distance;



	sys.stderr.write('\n');
	sys.stderr.write('Counting qnames present in sam_file2 that are missing from sam_file1...\n');
	num_processed = 0;
	for qname in sam_hash2.keys():
		num_processed += 1;
		if ((num_processed % 1000) == 0):
			sys.stderr.write('\rProcessed %d alignments...' % num_processed);
		sam_hash2_list = sam_hash2[qname];
		if (len(sam_hash2_list) > 0):
			if (sam_hash2_list[0].IsMapped() == True):
				num_mapped_2 += 1;

				try:
					sam_hash1_list = sam_hash1[qname];
					if (sam_hash1_list[0].IsMapped() == False):
						not_in_sam_file1 += 1;
						qnames_not_in_sam_file1.append([sam_hash2_list[0].evalue, qname]);
				except:
					not_in_sam_file1 += 1;
					qnames_not_in_sam_file1.append([sam_hash2_list[0].evalue, qname]);

		# if (len(sam_hash2_list) > 0):
		# 	if (sam_hash2_list[0].IsMapped() == True):
		# 		num_mapped_2 += 1;
		# 	if (qname in sam_hash1.keys()):
		# 		pass;
		# 	else:
		# 		not_in_sam_file1 += 1;
		# 		qnames_not_in_sam_file1.append(qname);
	sys.stderr.write('\n');
	sys.stderr.write('\n');

	fp_out = None;
	fp_out_lt0bp = None;
	fp_out_gt5000bp = None;
	out_file = out_summary_prefix + '.csv';
	out_file_lt0bp = out_summary_prefix + '_lt0bp.csv';
	out_file_gt5000bp = out_summary_prefix + '_gt5000bp.csv';
	
	if (out_summary_prefix != ''):
		try:
			fp_out = open(out_file, 'w');
			fp_out_lt0bp = open(out_file_lt0bp, 'w');
			fp_out_gt5000bp = open(out_file_gt5000bp, 'w');
		except IOError:
			sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_file));
			return;
			# exit(1);

	summary_line = '';
	summary_line += 'SAM file 1: %s\n' % sam_file1;
	summary_line += 'SAM file 2: %s\n' % sam_file2;
	summary_line += 'Number of qnames not present in SAM file 1: %d\n' % (not_in_sam_file1);
	summary_line += 'Number of qnames not present in SAM file 2: %d\n' % (not_in_sam_file2);
	summary_line += 'Number of qnames mapped to different references: %d\n' % (num_different_reference);
	summary_line += 'Number of alignments of different orientation: %d\n' % (num_different_orientation);
	summary_line += 'Number of shared qnames: %d\n' % (len(shared_qnames.keys()));
	summary_line += 'Mapped in SAM 1: %d\n' % (num_mapped_1);
	summary_line += 'Unmapped in SAM 1: %d\n' % (num_not_mapped_1);
	summary_line += 'Mapped in SAM 2: %d\n' % (num_mapped_2);
	summary_line += 'Unmapped in SAM 2: %d\n' % (num_not_mapped_2);
	summary_line += '\n';

	length_threshold = 9000;

	sys.stdout.write(summary_line);
	if (out_summary_prefix != ''):
		fp_out.write(summary_line);
	summary_line = '';

	summary_line_lt0bp = '';
	summary_line_gt5000bp = '';

	num_same_alignments = 0;
	i = 0;
	# while i < len(distance_to_qname_hash.keys()
	# print distance_to_qname_hash;
	for distance in sorted(distance_to_qname_hash.keys()):
		sorted_by_length = sorted(distance_to_sam_hash[distance], reverse=True, key=lambda sam_line: len(sam_line.seq));
		# sorted_qnames = ['%s <%d, %d>' % (single_sam_line.qname, len(single_sam_line.seq), single_sam_line.mapq) for single_sam_line in sorted_by_length];
		# positions = qname_to_pos[distance_to_qname_hash[distance]];
		# sorted_qnames = ['%s <len:%d, SAM1:%d, SAM2:%d>' % (single_sam_line.qname, len(single_sam_line.seq), positions[0], positions[1]) for single_sam_line in sorted_by_length];
		sorted_qnames = ['%s <len:%d, pos:%d>' % (single_sam_line.qname, len(single_sam_line.seq), single_sam_line.clipped_pos) for single_sam_line in sorted_by_length];

		sorted_qnames_above_length = [('%s' % (single_sam_line.qname)) for single_sam_line in sorted_by_length if (len(single_sam_line.seq) > length_threshold)];
		if (distance == 0):
			summary_line_lt0bp = ' \\\n'.join(sorted_qnames_above_length);
		if (distance > 5000):
			if (len(summary_line_gt5000bp) > 0):
				summary_line_gt5000bp += ' \\\n';
			summary_line_gt5000bp += ' \\\n'.join(sorted_qnames_above_length);

		# sorted_qnames = [str(len(single_sam_line.seq)) for single_sam_line in sorted(distance_to_sam_hash[distance], reverse=True, key=lambda sam_line: len(sam_line.seq))];
		# summary_line = str(distance) + '\t' + str(len(distance_to_qname_hash[distance])) + '\t' + '\t'.join(distance_to_qname_hash[distance]) + '\n';
		summary_line = str(distance) + '\t' + str(len(distance_to_qname_hash[distance])) + '\t' + '\t'.join(sorted_qnames) + '\n';
		if (distance <= distance_threshold):
			num_same_alignments += len(distance_to_qname_hash[distance]);

		# sys.stdout.write(summary_line);
		if (out_summary_prefix != ''):
			fp_out.write(summary_line);
		summary_line = '';

	summary_line = '\n';
	summary_line += 'Qnames that map to different references (total: %d):\n' % (len(qnames_different_ref));
	for qname_different_ref in qnames_different_ref:
		summary_line += '\t' + qname_different_ref + '\n';
	summary_line += '\n';
	summary_line += 'Qnames that map to different orientations (total: %d):\n' % (len(qnames_different_orient));
	for qname_different_orient in qnames_different_orient:
		summary_line += '\t' + qname_different_orient + '\n';
	summary_line += '\n';
	if (out_summary_prefix != ''):
		fp_out.write(summary_line);

	summary_line = 'Distance threshold to consider mappings same: %d\n' % distance_threshold;
	summary_line += 'Number of same mappings: %d\n' % num_same_alignments;
	summary_line += '(verbose) Number of same mappings: %d (%.2f%% in SAM1 / %.2f%% in SAM2) within %d bp distance.\n' % (num_same_alignments, 100.0 * float(num_same_alignments) / float(num_mapped_1 + num_not_mapped_1), 100.0 * float(num_same_alignments) / float(num_mapped_2 + num_not_mapped_2), distance_threshold);
	summary_line += '\n';
	sys.stdout.write(summary_line);
	if (out_summary_prefix != ''):
		fp_out.write(summary_line);
		fp_out_lt0bp.write(summary_line_lt0bp);
		fp_out_gt5000bp.write(summary_line_gt5000bp);
		summary_line = '';
		summary_line_lt0bp = '';
		summary_line_gt5000bp = '';

		sam1_basename = os.path.splitext(os.path.basename(sam_file1))[0];
		sam2_basename = os.path.splitext(os.path.basename(sam_file2))[0];

		out_file_qnames_only_in_sam2 = out_summary_prefix + '_qnames_only_in_%s.csv' % (sam2_basename);
		out_file_qnames_only_in_sam1 = out_summary_prefix + '_qnames_only_in_%s.csv' % (sam1_basename);
		out_file_qnames_only_in_sam2_as_sam = out_summary_prefix + '_qnames_only_in_%s.sam' % (sam2_basename);
		out_file_qnames_only_in_sam1_as_sam = out_summary_prefix + '_qnames_only_in_%s.sam' % (sam1_basename);
		out_file_qnames_in_both_sam1_as_sam = out_summary_prefix + '_qnames_in_both-alignments_from_%s.sam' % (sam1_basename);
		out_file_qnames_in_both_sam2_as_sam = out_summary_prefix + '_qnames_in_both-alignments_from_%s.sam' % (sam2_basename);

		summary_line += 'Output files:\n';
		summary_line += '\t%s\n' % (out_file_qnames_only_in_sam1);
		summary_line += '\t%s\n' % (out_file_qnames_only_in_sam2);
		summary_line += '\t%s\n' % (out_file_qnames_only_in_sam1_as_sam);
		summary_line += '\t%s\n' % (out_file_qnames_only_in_sam2_as_sam);
		summary_line += '\t%s\n' % (out_file_qnames_in_both_sam1_as_sam);
		summary_line += '\t%s\n' % (out_file_qnames_in_both_sam2_as_sam);

		try:
			fp_out_qnames_only_in_sam2 = open(out_file_qnames_only_in_sam2, 'w');
			fp_out_qnames_only_in_sam1 = open(out_file_qnames_only_in_sam1, 'w');
			# fp_out_qnames_only_in_sam2.write('\n'.join(qnames_not_in_sam_file1) + '\n');
			fp_out_qnames_only_in_sam2.write('\n'.join(['%e\t%s' % (value[0], value[1]) for value in sorted(qnames_not_in_sam_file1, key=lambda x: x[0])]) + '\n');
			fp_out_qnames_only_in_sam1.write('\n'.join(['%e\t%s' % (value[0], value[1]) for value in sorted(qnames_not_in_sam_file2, key=lambda x: x[0])]) + '\n');
			fp_out_qnames_only_in_sam2.close();
			fp_out_qnames_only_in_sam1.close();

			fp_out1 = open(out_file_qnames_only_in_sam2_as_sam, 'w');
			fp_out1.write('\n'.join(sam_headers2) + '\n');
			for value in sorted(qnames_not_in_sam_file1, key=lambda x: x[0]):
				fp_out1.write('\n'.join([sam_line.original_line for sam_line in sam_hash2[value[1]]]) + '\n');
			fp_out1.close();

			fp_out2 = open(out_file_qnames_only_in_sam1_as_sam, 'w');
			fp_out2.write('\n'.join(sam_headers1) + '\n');
			for value in sorted(qnames_not_in_sam_file2, key=lambda x: x[0]):
				fp_out2.write('\n'.join([sam_line.original_line for sam_line in sam_hash1[value[1]]]) + '\n');
			fp_out2.close();

			fp_out1 = open(out_file_qnames_in_both_sam1_as_sam, 'w');
			fp_out1.write('\n'.join(sam_headers1) + '\n');
			for value in shared_qnames:
				fp_out1.write('\n'.join([sam_line.original_line for sam_line in sam_hash1[value]]) + '\n');
			fp_out1.close();

			fp_out2 = open(out_file_qnames_in_both_sam2_as_sam, 'w');
			fp_out2.write('\n'.join(sam_headers2) + '\n');
			for value in shared_qnames:
				fp_out2.write('\n'.join([sam_line.original_line for sam_line in sam_hash2[value]]) + '\n');
			fp_out2.close();

		except IOError:
			sys.stderr.write('ERROR: Could not open file(s) for writing! Either "%s" or "%s".\n' % (out_file_qnames_only_in_sam2, out_file_qnames_only_in_sam1));

	if (out_summary_prefix != ''):
		fp_out.close();
		fp_out_lt0bp.close();
		fp_out_gt5000bp.close();



	# sys.stderr.write('Starting to count the number of correctly mapped bases in the tested SAM file!\n');

# src/sam_compare_cigars.py /home/ivan/work/eclipse-workspace/golden-bundle/alignments_for_testing/reads-simulated/PacBio-100k/escherichia_coli/graphmap-params_SSW_r-test-drive.sam /home/ivan/work/eclipse-workspace/golden-bundle/reads-simulated/PacBio-100k/escherichia_coli/reads.sam temp/cigcompare.temp
if __name__ == "__main__":
	if (len(sys.argv) < 4 or len(sys.argv) > 5):
		sys.stderr.write('Compares alignment positions between two SAM files.');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <input_sam_file_1> <input_sam_file_2> <out_file_prefix> [distance_threshold]\n' % sys.argv[0]);
		sys.stderr.write('\n');
		sys.stderr.write('\tdistance_threshold - default value is 100\n');
		sys.stderr.write('\n');
		sys.stderr.write('\tSeveral files will be created on disk:\n');
		sys.stderr.write('\t\t- <out_file_prefix>_qnames_only_in_sam2.csv\n');
		sys.stderr.write('\t\t- <out_file_prefix>_qnames_only_in_sam2.sam\n');
		sys.stderr.write('\t\t- <out_file_prefix>_qnames_in_both_sam1.sam\t- Alignments from SAM1 which are also mapped in SAM2.\n');
		sys.stderr.write('\t\t- <out_file_prefix>_qnames_only_in_sam1.csv\n');
		sys.stderr.write('\t\t- <out_file_prefix>_qnames_only_in_sam1.sam\n');
		sys.stderr.write('\t\t- <out_file_prefix>_qnames_in_both_sam2.sam\t- Alignments from SAM2 which are also mapped in SAM1.\n');
		sys.stderr.write('\t\t- <out_file_prefix>_lt0bp.csv\n');
		sys.stderr.write('\t\t- <out_file_prefix>_gt5000bp.csv\n');
		sys.stderr.write('\t\t- <out_file_prefix>.csv\n');
		sys.stderr.write('\n');
		exit(1);
	
	sam_file1 = sys.argv[1];
	sam_file2 = sys.argv[2];
	distance_threshold = 100;
	out_summary_prefix = sys.argv[3];

	if (len(sys.argv) == 5):
		distance_threshold = int(sys.argv[4]);
	
	out_summary_prefix = os.path.abspath(out_summary_prefix);
	if (not os.path.exists(os.path.dirname(out_summary_prefix))):
		sys.stderr.write('Creating output folder: "%s".\n' % (os.path.dirname(out_summary_prefix)));
		os.makedirs(os.path.dirname(out_summary_prefix));

	CompareTwoSAMs(sam_file1, sam_file2, distance_threshold, out_summary_prefix);
