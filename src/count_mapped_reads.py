#! /usr/bin/env python

import os;
import sys;
import utility_sam;

def CountMappedReads(sam_file, out_stats_prefix='', num_reads_in_reference_sam=1):
	fp = None;
	
	try:
		fp = open(sam_file, 'r');
	except IOError:
		print '[%s] ERROR: Could not open file "%s" for reading!' % ('CountMappedReads', sam_file);
		exit(1);
	
	num_alignments = 0;
	num_mapped_alignments = 0;
	sam_reads = {};
	sam_mapped_reads = {};
	num_mapped_bases = 0;	# Count of the number of non-clipped bases in an alignment. Only highest-scoring alignments are considered for each read,
	total_num_bases = 0;
	highest_scoring_alignment_for_read = {};	# Hash to keep track of the highest scoring alignment for a read.
	
	for line in fp:
		if (len(line) == 0 or line[0] == '@'):
			continue;
		
		sam_line = utility_sam.SAMLine(line.rstrip());
		
		# Count the number of occurances of each read in the SAM file (i.e. multiple mappings).
		try:
			sam_reads[sam_line.qname] += 1;
			#print sam_line.qname;
		except:
			sam_reads[sam_line.qname] = 1;
			
		if (sam_line.IsMapped() == True):
			# Count the occurances of mapped reads.
			try:
				sam_mapped_reads[sam_line.qname] += 1;
			except:
				sam_mapped_reads[sam_line.qname] = 1;
			
			try:
				if (highest_scoring_alignment_for_read[sam_line.qname].chosen_quality < sam_line.chosen_quality):
					highest_scoring_alignment_for_read[sam_line.qname] = sam_line;
			except:
				highest_scoring_alignment_for_read[sam_line.qname] = sam_line;
			
			num_mapped_alignments += 1;
		
		num_alignments += 1;
	
	fp.close();
	
	num_unique_reads = 0;
	for value in sam_reads.values():
		if (value == 1):
			num_unique_reads += 1;
	num_unique_mapped_reads = 0;
	for value in sam_mapped_reads.values():
		if (value == 1):
			num_unique_mapped_reads += 1;
	
	for best_read in highest_scoring_alignment_for_read.values():
		total_num_bases += len(sam_line.seq);
		num_mapped_bases += best_read.CalcNumMappedBases();
	
	
	
	ret_lines = '';
	ret_lines += 'Output prefix: %s\n' % out_stats_prefix;
	ret_lines += 'Number of alignments in SAM file: %d\n' % num_alignments;
	ret_lines += 'Number of alignments in SAM file which are marked mapped: %d (%.2f%%)\n' % (num_mapped_alignments, (float(num_mapped_alignments) / float(num_alignments))*100.0);
	ret_lines += 'Number of reads in SAM file: %d (%.2f%%)\n' % (len(sam_reads.keys()), (float(len(sam_reads.keys())) / float(num_alignments))*100.0);
	ret_lines += 'Number of mapped reads in SAM file: %d (percent_reads: %.2f%%, percent_alignments: %.2f%%)\n' % (len(sam_mapped_reads.keys()), (float(len(sam_mapped_reads.keys()) / float(len(sam_reads.keys())))*100.0), (float(len(sam_mapped_reads.keys())) / float(num_alignments))*100.0);
	ret_lines += 'Number of unique reads in SAM file: %d (percent_reads: %.2f%%, percent_alignments: %.2f%%)\n' % (num_unique_reads, (float(num_unique_reads) / float(len(sam_reads.keys())))*100.0, (float(num_unique_reads) / float(num_alignments))*100.0);
	ret_lines += 'Number of reads with multiple alignments: %d (percent_reads: %.2f%%, percent_alignments: %.2f%%)\n' % ((len(sam_reads.keys()) - num_unique_reads), (float((len(sam_reads.keys()) - num_unique_reads)) / float(len(sam_reads.keys())))*100.0, (float((len(sam_reads.keys()) - num_unique_reads)) / float(num_alignments))*100.0);
	ret_lines += 'Number of reads in reference SAM: %d\n' % (num_reads_in_reference_sam);
	ret_lines += 'Percent of reads compared to the reference SAM: %.2f%%\n' % ((float(len(sam_reads.keys())) / float(num_reads_in_reference_sam)) * 100.0);
	ret_lines += 'Percent of mapped reads compared to the reference SAM: %.2f%%\n' % ((float(len(sam_mapped_reads.keys())) / float(num_reads_in_reference_sam)) * 100.0);
	ret_lines += 'Percent of reads with unique alignments compared to the reference SAM: %.2f%%\n' % ((float(num_unique_reads) / float(num_reads_in_reference_sam)) * 100.0);
	ret_lines += 'Percent of reads with multiple alignments compared to the reference SAM: %.2f%%\n' % ((float(len(sam_reads.keys()) - num_unique_reads) / float(num_reads_in_reference_sam)) * 100.0);
	ret_lines += 'Number of mapped bases: %ld\n' % (num_mapped_bases);
	ret_lines += 'Total number of bases: %ld\n' % (total_num_bases);
	ret_lines += 'Percent of mapped bases: %.2f\n' % ((float(num_mapped_bases) / float(total_num_bases)) * 100.0);

	ret_plot_lines = '';
	ret_plot_lines += '%s\n' % sam_file;
	ret_plot_lines += '%s\n' % out_stats_prefix;
	ret_plot_lines += 'x\t';
	ret_plot_lines += '-\t';
	ret_plot_lines += 'counts\t';
	ret_plot_lines += 'num_alignments_in_sam\t';
	ret_plot_lines += 'num_alignments_in_sam_marked_mapped\t'
	ret_plot_lines += 'num_reads_in_sam\t';
	ret_plot_lines += 'num_reads_in_sam_marked_mapped\t';
	ret_plot_lines += 'num_reads_with_unique_alignments_in_sam\t';
	ret_plot_lines += 'num_reads_with_multiple_alignments_in_sam\t';
	ret_plot_lines += 'num_reads_in_reference_sam\t';
	ret_plot_lines += 'percent_reads_in_sam_compared_to_reference_sam\t';
	ret_plot_lines += 'percent_mapped_reads_in_sam_compared_to_reference_sam\t';
	ret_plot_lines += 'percent_reads_with_unique_alignments_in_sam_compared_to_reference_sam\t';
	ret_plot_lines += 'percent_reads_with_multiple_alignments_in_sam_compared_to_reference_sam\t';
	ret_plot_lines += 'num_mapped_bases';
	ret_plot_lines += 'total_num_bases';
	ret_plot_lines += 'percent_mapped_bases';
	ret_plot_lines += '\n';
	
	ret_plot_lines += 'y\t';
	ret_plot_lines += '%s\t' % os.path.splitext(os.path.basename(sam_file))[0];
	ret_plot_lines += 'N\t';
	ret_plot_lines += '%d\t' % (num_alignments);
	ret_plot_lines += '%d\t' % (num_mapped_alignments);
	ret_plot_lines += '%d\t' % (len(sam_reads.keys()));
	ret_plot_lines += '%d\t' % (len(sam_mapped_reads.keys()));
	ret_plot_lines += '%d\t' % (num_unique_reads);
	ret_plot_lines += '%d\t' % (len(sam_reads.keys()) - num_unique_reads);
	ret_plot_lines += '%s\t' % (num_reads_in_reference_sam);
	ret_plot_lines += '%.2f\t' % ((float(len(sam_reads.keys())) / float(num_reads_in_reference_sam)) * 100.0);
	ret_plot_lines += '%.2f\t' % ((float(len(sam_mapped_reads.keys())) / float(num_reads_in_reference_sam)) * 100.0);
	ret_plot_lines += '%.2f\t' % ((float(num_unique_reads) / float(num_reads_in_reference_sam)) * 100.0);
	ret_plot_lines += '%.2f\t' % ((float(len(sam_reads.keys()) - num_unique_reads) / float(num_reads_in_reference_sam)) * 100.0);
	ret_plot_lines += '%d\t' % (num_mapped_bases);
	ret_plot_lines += '%ld\t' % (total_num_bases);
	ret_plot_lines += '%.2f' % ((float(num_mapped_bases) / float(total_num_bases)) * 100.0);
	ret_plot_lines += '\n';
	
	WriteSAMStats(ret_lines, ret_plot_lines, out_stats_prefix);
	
	#ret_plot_lines += 'num_alignments_in_sam\t%d\n' % num_alignments;
	#ret_plot_lines += 'num_alignments_in_sam_marked_mapped\t%d\n' % num_mapped_alignments;
	#ret_plot_lines += 'num_reads_in_sam\t%d\n' % (len(sam_reads.keys()));
	#ret_plot_lines += 'num_mapped_reads_in_sam\t%d\n' % (len(sam_mapped_reads.keys()));
	#ret_plot_lines += 'num_reads_with_unique_alignments_in_sam\t%d\n' % (num_unique_reads);
	#ret_plot_lines += 'num_reads_with_multiple_alignments_in_sam\t%d\n' % (len(sam_reads.keys()) - num_unique_reads);


	#print ret_lines;
	
	return [ret_lines, ret_plot_lines];

def WriteSAMStats(readable_lines, plot_lines, out_stats_prefix):
	if (out_stats_prefix != ''):
		try:
			fp = open(out_stats_prefix + '.txt', 'w');
			fp.write(readable_lines);
			fp.close();
		except Exception, e:
			print 'ERROR: Could not open file "%s" for writing!' % (out_stats_prefix + '.txt');
			print e;
			
		try:
			fp = open(out_stats_prefix + '.plot', 'w');
			fp.write(plot_lines);
			fp.close();
		except Exception, e:
			print 'ERROR: Could not open file "%s" for writing!' % (out_stats_prefix + '.plot');
			print e;
	
def CollectSAMStats(stats_prefix, suppress_error_messages=False):
	readable_lines = '';
	plot_lines = '';
	
	try:
		fp = open(stats_prefix + '.txt', 'r');
		lines = fp.readlines();
		fp.close();
		
		readable_lines = ''.join(lines);
	except IOError:
		if (suppress_error_messages == False):
			print 'ERROR: Could not open file "%s" for reading!' % (stats_prefix + '.txt');
	
	try:
		fp = open(stats_prefix + '.plot', 'r');
		lines = fp.readlines();
		fp.close();
		
		plot_lines = ''.join(lines);
	except IOError:
		if (suppress_error_messages == False):
			print 'ERROR: Could not open file "%s" for reading!' % (stats_prefix + '.plot');
	
	return [readable_lines, plot_lines];

if __name__ == "__main__":
	
	if (len(sys.argv) < 2):
		print 'Tool for counting the number of reads and the number of mapped bases in a given SAM file.'
		print 'Every mapped read is counted only once (even though it might have multiple alignments).'
		print 'In case of mapped-base counting, only one alignment (the one with highest MAPQ/AS score)';
		print 'is chosen.';
		print ' ';
		print 'Usage:';
		print '\t%s <sam_file> [<output_stats_prefix>] [num_reads_in_reference_sam]' % sys.argv[0];
		exit(1);
	
	sam_file = sys.argv[1];
	out_stats_prefix = '';
	num_reads_in_reference_sam = 1;
	if (len(sys.argv) >= 3):
		out_stats_prefix = sys.argv[2];
	if (len(sys.argv) >= 4):
		num_reads_in_reference_sam = int(sys.argv[3]);
	[ret_lines, ret_plot_lines] = CountMappedReads(sam_file, out_stats_prefix, num_reads_in_reference_sam);
	print ret_lines;
