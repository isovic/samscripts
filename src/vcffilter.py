#! /usr/bin/env python

import re;
import os;
import sys;

def extract_region(vcf_file, chrname, start_pos, end_pos):
	fp_in = None;
	
	try:
		fp_in = open(vcf_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, vcf_file));
		exit(1);
	
	num_accepted = 0;
	num_rejected = 0;

	fp_out = sys.stdout;

	pos_position = -1;
	
	i = 0;
	for line in fp_in:
		i += 1;
		if ((i % 1000) == 0):
			sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		if (len(line.strip()) == 0 or line[0] == '#'):
			fp_out.write(line);

			if (len(line) > 1 and line[1] != '#'):
				params = line[1:].split('\t');
				current_param = 0;
				for param in params:
					if (param.lower() == 'pos'):
						pos_position = current_param;
						break;
					current_param += 1;
			continue;

		split_line = line.strip().split('\t');
		
		if ((i % 1000) == 0):
			sys.stderr.write(', chr: "%s"' % (split_line[0]));

		if (split_line[0] == chrname):
			if (pos_position < 0):
				sys.stderr.write('ERROR: Could not find parameter description line, needed to determine the POS parameter!\n');
				exit(1);

			pos = int(split_line[pos_position]);

			if (pos >= start_pos and pos <= end_pos):
				fp_out.write(line);
				num_accepted += 1;
			else:
				num_rejected += 1;
		else:
			num_rejected += 1;

	fp_in.close();
	
	sys.stderr.write('\n');
	sys.stderr.write('num_accepted = %d (%.2f%%)\n' % (num_accepted, (float(num_accepted) / float(num_accepted + num_rejected)) * 100.0));
	sys.stderr.write('num_rejected = %d (%.2f%%)\n' % (num_rejected, (float(num_rejected) / float(num_accepted + num_rejected)) * 100.0));

def split_multiallelic_snps(vcf_file):
	fp_in = None;
	
	try:
		fp_in = open(vcf_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, vcf_file));
		exit(1);

	fp_out = sys.stdout;

	alt_position = -1;
	
	i = 0;
	for line in fp_in:
		i += 1;
		if ((i % 1000) == 0):
			sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		if (len(line.strip()) == 0 or line[0] == '#'):
			fp_out.write(line);

			if (len(line) > 1 and line[1] != '#'):
				params = line[1:].split('\t');
				current_param = 0;
				for param in params:
					if (param.lower() == 'alt'):
						alt_position = current_param;
						break;
					current_param += 1;
			continue;

		split_line = line.strip().split('\t');

		if ((i % 1000) == 0):
			sys.stderr.write(', chr: "%s"' % (split_line[0]));

		if (alt_position < 0):
			sys.stderr.write('ERROR: Could not find parameter description line, needed to determine the ALT parameter!\n');
			exit(1);

		if ((',' in split_line[alt_position]) == False):
			fp_out.write(line);
		else:
			split_alt = split_line[alt_position].split(',');
			for alt in split_alt:
				split_line[alt_position] = alt;
				fp_out.write('\t'.join(split_line) + '\n');

	fp_in.close();
	
	sys.stderr.write('\n');

def vcf_sort(vcf_file):
	fp_in = None;
	
	try:
		fp_in = open(vcf_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, vcf_file));
		exit(1);

	fp_out = sys.stdout;

	pos_position = -1;
	
	vcf_header = [];
	vcf_lines = [];

	i = 0;
	for line in fp_in:
		i += 1;
		if ((i % 1000) == 0):
			sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		if (len(line.strip()) == 0 or line[0] == '#'):
			# fp_out.write(line);
			vcf_header.append(line.strip());

			if (len(line) > 1 and line[1] != '#'):
				params = line[1:].split('\t');
				current_param = 0;
				for param in params:
					if (param.lower() == 'pos'):
						pos_position = current_param;
						break;
					current_param += 1;
			continue;

		split_line = line.strip().split('\t');

		if ((i % 1000) == 0):
			sys.stderr.write(', chr: "%s"' % (split_line[0]));

		if (pos_position < 0):
			sys.stderr.write('ERROR: Could not find parameter description line, needed to determine the POS parameter!\n');
			exit(1);

		vcf_lines.append([int(split_line[pos_position]), line.strip()]);

	fp_in.close();
	
	sys.stderr.write('\n');

	sorted_vcf_lines = sorted(vcf_lines, key=lambda x: x[0]);
	fp_out.write('\n'.join(vcf_header) + '\n');
	fp_out.write('\n'.join([line[1] for line in sorted_vcf_lines]) + '\n');

def count_nonpass_variants(vcf_file, verbose=True):
	fp_in = None;
	
	try:
		fp_in = open(vcf_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!' % (__name__, vcf_file));
		exit(1);

	fp_out = sys.stdout;

	filter_position = -1;
	info_position = -1;

	num_pass_snps = 0;
	num_nonpass_snps = 0;

	i = 0;
	for line in fp_in:
		i += 1;
		if ((i % 1000) == 0):
			sys.stderr.write('\rLine %d, num_accepted: %d, num_rejected: %d' % (i, num_accepted, num_rejected));

		line = line.strip();
		if (len(line) == 0 or line[0] == '#'):

			if (len(line) > 1 and line[1] != '#'):
				params = line[1:].split('\t');
				current_param = 0;
				for param in params:
					if (param.lower() == 'filter'):
						filter_position = current_param;
					if (param.lower() == 'info'):
						info_position = current_param;
					current_param += 1;
			continue;

		split_line = line.strip().split('\t');

		if ((i % 1000) == 0):
			sys.stderr.write(', chr: "%s"' % (split_line[0]));

		if (filter_position < 0 or info_position < 0):
			sys.stderr.write('ERROR: Could not find parameter description line, needed to determine the POS parameter!\n');
			exit(1);

		# vcf_lines.append([int(split_line[pos_position]), line.strip()]);
		split_info = split_line[info_position].split(';');
		info_vals = {};
		for value in split_info:
			split_value = value.split('=');
			if (len(split_value) != 2):
				continue;
			try:
				info_vals[split_value[0].lower()] = split_value[1].lower();
			except Exception, e:
				sys.stderr.write('value = "%s"\n' % value);
				sys.stderr.write(str(e));
				exit(1);
		if ('vartype' in info_vals.keys()):
			vartype = info_vals['vartype'];
		elif ('type' in info_vals.keys()):
			vartype = info_vals['type'];
		else:
			# sys.stderr.write('ERROR: VCF line does not contain varType info! Counting it as a SNP.\n');
			# sys.stderr.write(line);
			# sys.stderr.write(str(info_vals.keys()) + '\n');
			vartype = 'snp';
			# continue;

		if (vartype.lower().split(',')[0] != 'snp'):
			continue;

		if (split_line[filter_position].lower() == 'pass'):
			# print '[%d, %d] PASS: %s' % (num_pass_snps, num_nonpass_snps, line.strip());
			num_pass_snps += 1;
		else:
			# print '[%d, %d] Non-PASS: %s' % (num_pass_snps, num_nonpass_snps, line.strip());
			num_nonpass_snps += 1;

	fp_in.close();

	if (verbose == True):
		sys.stdout.write('%d\t%d\n' % (num_pass_snps, num_nonpass_snps));

	return [num_pass_snps, num_nonpass_snps];



if __name__ == "__main__":
	if (len(sys.argv) < 2):
		sys.stderr.write('Various filtering methods for filtering VCF files.\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\tregion\n');
		sys.stderr.write('\tsplitsnps\n');
		sys.stderr.write('\tsort\n');
		sys.stderr.write('\tcountuncertsnps\n');

		exit(0);

	if (sys.argv[1] == 'region'):
		if (len(sys.argv) != 6):
			sys.stderr.write('Extracts only variant lines which fall within a defined region. Output is on stdout. Comment lines are output as well.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_vcf_file> chrname startpos endpos\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		vcf_file = sys.argv[2];
		chrname = sys.argv[3];
		start_pos = int(sys.argv[4]);
		end_pos = int(sys.argv[5]);
		extract_region(vcf_file, chrname, start_pos, end_pos);
		exit(0);

	elif (sys.argv[1] == 'splitsnps'):
		if (len(sys.argv) != 3):
			sys.stderr.write('Splits multiallelic SNP variants into separate lines.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_vcf_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		vcf_file = sys.argv[2];
		split_multiallelic_snps(vcf_file);
		exit(0);

	elif (sys.argv[1] == 'sort'):
		if (len(sys.argv) != 3):
			sys.stderr.write('Sorts a VCF file by position.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_vcf_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		vcf_file = sys.argv[2];
		vcf_sort(vcf_file);
		exit(0);

	elif (sys.argv[1] == 'countuncertsnps'):
		if (len(sys.argv) != 3):
			sys.stderr.write('Count SNP variants with non PASS filter value.\n');
			sys.stderr.write('Usage:\n');
			sys.stderr.write('\t%s %s <input_vcf_file>\n' % (sys.argv[0], sys.argv[1]));
			exit(0);

		vcf_file = sys.argv[2];
		count_nonpass_variants(vcf_file);
		exit(0);

	else:
		sys.stderr.write('ERROR: Unknown subcommand!\n');
		exit(0);
