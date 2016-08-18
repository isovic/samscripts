#! /usr/bin/env python

# Copyright Ivan Sovic, 2015. www.sovic.org
#
# A script that takes an input FASTA/FASTQ file and performs various methods
# of filtering on it. Examples include: selecting sequences that contain certain
# strings in their header, and others.

import re
import os
import sys
import fastqparser
import numpy as np
import random
import operator;



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
                break;

    sys.stderr.write('\n');
    fp_in.close();

def filter_seqs_by_read_id(input_fastq_path, read_id_path, out_fastq_path, fp_out):
    try:
        fp_in = open(input_fastq_path, 'r');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
        exit(0);

    try:
        fp_filter = open(read_id_path, 'r');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % read_id_path);
        exit(0);

    filter_read_ids = fp_filter.readlines();
    fp_filter.close();
    filter_read_ids = [int(line.strip()) for line in filter_read_ids if (len(line.strip()) > 0)];
    id_hash = {};
    for read_id in filter_read_ids:
        id_hash[read_id] = 1;

    num_matches = 0;

    num_reads = 0;
    while True:
        [header, read] = fastqparser.get_single_read(fp_in);

        if (len(read) == 0):
            break;
        if (num_reads in id_hash):
            num_matches += 1;
            sys.stderr.write('\rFound %d seqs.' % (num_matches));
            fp_out.write('\n'.join(read) + '\n');

        num_reads += 1;

    sys.stderr.write('\n');
    fp_in.close();

def filter_seqs_by_header_list(input_fastq_path, filter_headers, out_fastq_path, fp_out):
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

        for filter_header in filter_headers:
            if ((filter_header[-1] == ';' and filter_header[0:-1].lower() == header.lower()) or
                (filter_header[-1] != ';' and filter_header.lower() in header.lower())):
                num_matches += 1;
                sys.stderr.write('\rFound %d seqs, last: "%s".' % (num_matches, header));
                fp_out.write('\n'.join(read) + '\n');
                break;

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

    i = 0;
    while True:
        i += 1;
        [header, read] = fastqparser.get_single_read(fp_in);

        if (len(read) == 0):
            break;

        header_id = header.split()[0].lower();

        if (header_id in id_hash):
            continue;

        num_matches += 1;
        if ((i % 1000) == 0):
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

def filter_seqs_by_length(input_fastq_path, length_threshold, is_less_than, out_fastq_path, fp_out):
    try:
        fp_in = open(input_fastq_path, 'r');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
        exit(0);

    num_matches = 0;

    i = 0;
    while True:
        i += 1;
        [header, read] = fastqparser.get_single_read(fp_in);

        if (len(read) == 0):
            break;

        if ((is_less_than == False and len(read[1]) >= length_threshold) or (is_less_than == True and len(read[1]) <= length_threshold)):
            num_matches += 1;
            if ((i % 1000) == 0):
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
        if ((num_reads % 100) == 0):
            sys.stderr.write('\rProcessing seq %d...' % num_reads);

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
        else:
            means_1d.append(np.mean(phreds));

        num_reads += 1;

    sys.stderr.write('\n');

    if (len(means_1d) > 0):
        sys.stdout.write('[1d] avg = %5.2f\tstd = %5.2f\tmin = %2d\tmax = %2d\n' % (np.mean(means_1d), np.std(means_1d), min(means_1d), max(means_1d)));
    else:
        sys.stdout.write('[1d] avg = %5.2f\tstd = %5.2f\tmin = %2d\tmax = %2d\n' % (0.0, 0.0, 0.0, 0.0));
    if (len(means_2d) > 0):
        sys.stdout.write('[2d] avg = %5.2f\tstd = %5.2f\tmin = %2d\tmax = %2d\n' % (np.mean(means_2d), np.std(means_2d), min(means_2d), max(means_2d)));
    else:
        sys.stdout.write('[2d] avg = %5.2f\tstd = %5.2f\tmin = %2d\tmax = %2d\n' % (0.0, 0.0, 0.0, 0.0));

    fp_in.close();

def base_quality_filter(input_fastq_path, lte_gte, qv_threshold, out_fastq_path, fp_out):
    print 'lte_gte = "%s"' % (lte_gte);

    try:
        fp_in = open(input_fastq_path, 'r');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
        return;

    num_reads = 0;
    num_outputted_reads = 0;
    num_skipped_reads = 0;

    while True:
        if ((num_reads % 1000) == 0):
            sys.stderr.write('\rProcessing seq %d, (%d passed)...' % (num_reads, num_outputted_reads));

        [header, read] = fastqparser.get_single_read(fp_in);
        if (len(read) == 0):
            break;

        if (len(read) < 4):
            sys.stderr.write('Given file is not a FASTQ file! Exiting.\n');
            return;

        quals = read[3];
        phreds = [];
        for char in quals:
            phreds.append(ord(char) - 33);
        mean_qv = np.mean(phreds);

        if ((lte_gte == 'gt' and mean_qv > qv_threshold) or
            (lte_gte == 'gte' and mean_qv >= qv_threshold) or
            (lte_gte == 'lt' and mean_qv < qv_threshold) or
            (lte_gte == 'lte' and mean_qv <= qv_threshold) or
            (lte_gte == 'eq' and mean_qv == qv_threshold)):

            fp_out.write('\n'.join(read) + '\n');
            num_outputted_reads += 1;
        else:
            num_skipped_reads += 1;

        num_reads += 1;

    sys.stderr.write('\n');
    sys.stderr.write('Total number of sequences: %d\n' % (num_reads));
    sys.stderr.write('Number of outputted sequences: %d\n' % (num_outputted_reads));
    sys.stderr.write('Number of sequences not satisfying condition: %d\n' % (num_skipped_reads));
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
        if ((num_reads % 1000) == 0):
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

        read[1] = re.sub("(.{%d})"%(wrap_length), "\\1\n", read[1], 0, re.DOTALL);    ### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
        if (len(read) == 4):
            read[3] = re.sub("(.{%d})"%(wrap_length), "\\1\n", read[3], 0, re.DOTALL);    ### Wrap the qual line, because DALIGNER has a 9998bp line len limit.
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
        if (len(read) == 4):
            read[2] = '+' + str(current_read);
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
    # return header_hash;

def count_1d2d(input_fastq_path):
    try:
        fp_in = open(input_fastq_path, 'r');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
        exit(0);

    num_reads = 0;
    num_1d = 0;
    num_2d = 0;

    while True:
        if ((num_reads % 100) == 0):
            sys.stderr.write('\rProcessing seq %d...' % num_reads);

        [header, read] = fastqparser.get_single_read(fp_in);
        if (len(read) == 0):
            break;

        if (('twodir' in header.lower()) or ('2d' in header.lower())):
            num_2d += 1;
        else:
            num_1d += 1;

        num_reads += 1;

    sys.stderr.write('\n');

    fp_in.close();

    try:
        sys.stdout.write('num_1d = %d (%.2f%%)\n' % (num_1d, (float(num_1d) / float(num_reads) * 100.0)));
        sys.stdout.write('num_2d = %d (%.2f%%)\n' % (num_2d, (float(num_2d) / float(num_reads) * 100.0)));
        sys.stdout.write('num_reads = %d\n' % (num_reads));
    except Exception, e:
        sys.stderr.write(str(e));
        pass;



def check_nanopore_paths(input_fastq_path, fast5_root_path):
    try:
        fp_in = open(input_fastq_path, 'r');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
        exit(0);

    num_reads = 0;
    num_wrong_paths = 0;

    while True:
        if ((num_reads % 1000) == 0):
            sys.stderr.write('\rProcessing seq %d...' % num_reads);

        [header, read] = fastqparser.get_single_read(fp_in);
        if (len(read) == 0):
            break;

        path_from_header = header.split()[-1];

        if (os.path.exists('%s/%s' % (fast5_root_path, path_from_header)) == False):
            num_wrong_paths += 1;
            sys.stderr.write('[%d] ' % (num_wrong_paths));
            sys.stdout.write('%s\n' % (header));

        num_reads += 1;

    sys.stderr.write('\n');

    fp_in.close();

### Subsample
# Subsampling calculates a total number of bases that should be in the result (according to a given coverage)
# and then current total number of bases
# Uses the ratio of the two to randomly select which row will be output as a result.
# Depending on the randomly generated numbers, the exact coverage of the result could vary,
# but should be close enough to a give parameter, given a large enough sample
def subsample(input_fastq_path, desired_coverage, ref_genome_size):

    if input_fastq_path.endswith('.fa') or input_fastq_path.endswith('.fasta'):
        type = 'fasta'
    elif input_fastq_path.endswith('.fq') or input_fastq_path.endswith('.fastq'):
        type = 'fastq'
    else:
        sys.stderr.write('\nUnrecognized file extension! Assuming fasta format!')
        type = 'fasta'

    [headers, seqs, quals] = fastqparser.read_fastq(input_fastq_path)

    numBases = sum([len(seq) for seq in seqs])
    desiredNumBases = desired_coverage * ref_genome_size

    ratio = float(desiredNumBases) / numBases

    if desiredNumBases > numBases:
        sys.stderr.write('Insufficient data for coverage %d\n' % desired_coverage)
        exit(1)

    sys.stderr.write('\n\nStatistics:')
    sys.stderr.write('\nNumber of desired/current bases: %i/%i' % (desiredNumBases, numBases))
    sys.stderr.write('\nGenome size: %i' % ref_genome_size)
    sys.stderr.write('\nDesired coverage: %i' % desired_coverage)
    sys.stderr.write('\nRatio: %f' % ratio)

    numOutputBases = 0
    numOutputSeqs = 0
    for i in xrange(len(headers)):
        rnd = random.random()
        if rnd <= ratio:
            # Print the correspodning row
            numOutputBases += len(seqs[i])
            numOutputSeqs += 1
            if type == 'fasta':
                sys.stdout.write('>%s\n%s\n' % (headers[i], seqs[i]))
            elif type == 'fastq':
                sys.stdout.write('@%s\n%s\n+%s\n%s\n' % (headers[i], seqs[i], headers[i], quals[i]))
            else:
                exception('Critical Error! Invalid filetype')

    sys.stderr.write('\nOutput bases: %i' % numOutputBases)
    sys.stderr.write('\nOutput seqs: %i' % numOutputSeqs)
    sys.stderr.write('\n')


def getPaired(input_fastq_path, target_fastq_path):

    if input_fastq_path.endswith('.fa') or input_fastq_path.endswith('.fasta'):
        type = 'fasta'
    elif input_fastq_path.endswith('.fq') or input_fastq_path.endswith('.fastq'):
        type = 'fastq'
    else:
        sys.stderr.write('\nUnrecognized file extension! Assuming fasta format!')
        type = 'fasta'

    [iheaders, iseqs, iquals] = fastqparser.read_fastq(input_fastq_path)
    [theaders, tseqs, tquals] = fastqparser.read_fastq(target_fastq_path)
    sys.stderr.write('\n\nLoaded input and target files!\n')
    sys.stderr.write('ilen = %d, tlen = %d \n' % (len(iheaders), len(theaders)))

    sys.stderr.write('Hashing target file headers ...')
    tdict = {}
    for theader in theaders:
        theader = theader.split()[0]
        tdict[theader] = 1

    for i in xrange(len(iheaders)):
        iheader = iheaders[i]
        iheader = iheader.split()[0]
        if iheader in tdict:
            if type == 'fasta':
                sys.stdout.write('>%s\n%s\n' % (iheaders[i], iseqs[i]))
            elif type == 'fastq':
                sys.stdout.write('@%s\n%s\n+%s\n%s\n' % (iheaders[i], iseqs[i], iheaders[i], iquals[i]))
            else:
                exception('Critical Error! Invalid filetype')
#        else:
#            import pdb
#            pdb.set_trace()

        if (i%100 == 0):
            sys.stderr.write('\nProcessed %d reads' % i)


def fastq2fasta(input_fastq_path, replace_semicolon=True):

    if not (input_fastq_path.endswith('.fq') or input_fastq_path.endswith('.fastq')):
        sys.stderr.write('\nUnrecognized file extension! Assuming fastq or fasta format!')

    [headers, seqs, quals] = fastqparser.read_fastq(input_fastq_path)

    for i in xrange(len(headers)):
        if (replace_semicolon == True):
            newheader = headers[i].replace(':', ' ')
            sys.stdout.write('>%s\n%s\n' % (newheader, seqs[i]))
        else:
            sys.stdout.write('>%s\n%s\n' % (headers[i], seqs[i]))

def length_distribution(input_path):
    if (input_path.endswith('.fq') or input_path.endswith('.fastq')):
        pass
    elif (input_path.endswith('.fa') or input_path.endswith('.fasta')):
        pass
    else:
        sys.stderr.write('\nUnrecognized file extension! Assuming fastq or fasta format!')

    [headers, seqs, quals] = fastqparser.read_fastq(input_path)

    distrib = {}
    for i in xrange(len(headers)):
        length = len(seqs[i])
        if length in distrib:
            distrib[length] += 1
        else:
            distrib[length] = 1

    sys.stdout.write('\n\n%d sequences!' % len(distrib))
    for length,count in distrib.iteritems():
        sys.stdout.write('\n%d sequences of length %d' % (count, length))
    sys.stdout.write('\n')

def extract_subseqs(input_fastq_path, start_coord, end_coord, out_fastq_path, fp_out):
    try:
        fp_in = open(input_fastq_path, 'r');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % input_fastq_path);
        exit(0);

    num_matches = 0;

    i = 0;
    while True:
        i += 1;
        [header, read] = fastqparser.get_single_read(fp_in);

        if (len(read) == 0):
            break;

        if (start_coord < len(read[1])):
            current_end = end_coord if (end_coord < len(read[1])) else len(read[1]);
            read[1] = read[1][start_coord:current_end];
            if (len(read) == 4):
                read[3] = read[3][start_coord:current_end];
            fp_out.write('\n'.join(read) + '\n');

        num_matches += 1;
        if ((i % 1000) == 0):
            sys.stderr.write('\rProcessed %d seqs.' % (num_matches));

    sys.stderr.write('\n');
    fp_in.close();

def msa2fasta(input_path, fp_out):
    [headers, seqs, quals] = fastqparser.read_fastq(input_path)
    for i in xrange(0, len(seqs)):
        seqs[i] = seqs[i].upper();

    cons_seq = '';
    for i in xrange(0, len(seqs[0])):
        base_counts = {'A': 0, 'C': 0, 'T': 0, 'G': 0, '.': 0, '-': 0};
        for j in xrange(0, len(seqs)):
            base_counts[seqs[j][i]] += 1;
        sorted_base_counts = sorted(base_counts.items(), key=operator.itemgetter(1));
        # Print sorted_base_counts;
        if ((sorted_base_counts[-1][0] in '.-') == False):
            cons_seq += sorted_base_counts[-1][0]
    fp_out.write('>Consensus_from_MSA\n%s\n' % (cons_seq));

def separate_seqs(input_fastq_file, out_folder):
    if (not os.path.exists(out_folder)):
        os.makedirs(out_folder);

    [headers, seqs, quals] = fastqparser.read_fastq(input_fastq_file);

    for i in xrange(len(seqs)):
        fp = open('%s/%d.fast%c' % (out_folder, (i+1), input_fastq_file[-1]), 'w');
        # fp = open('%s/%d' % (out_folder, (i+1)), 'w');
        if (len(quals[i]) == 0):
            fp.write('>%s\n%s\n' % (headers[i], seqs[i]));
        else:
            fp.write('@%s\n%s\n+\n%s\n' % (headers[i], seqs[i], quals[i]));
    fp.close();

def convert_reads_to_pacbio_format(reads_file, pacbio_reads_file):
    try:
        fp_in = open(reads_file, 'r');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % reads_file);
        exit(0);

    try:
        fp_out = open(pacbio_reads_file, 'w');
    except:
        sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % pacbio_reads_file);
        exit(0);

    current_read = 0;

    header_conversion_hash = {};

    while True:
        [header, read] = fastqparser.get_single_read(fp_in);
        
        if (len(read) == 0):
            break;

        current_read += 1;

        if (len(read[1]) <= 50):    ### DALIGNER has a lower length limit of 'k'. Trim more just to be safe.
            sys.stderr.write('Found a read shorter than 10bp. Removing from the output.\n');
            continue;

        ### Check if the read is already formatted like PacBio.
        if (header.count('/') == 2 and 'RQ' in header):
            fp_out.write('\n'.join(read) + '\n');
            continue;

        trimmed_header = header.replace('_', ' ').split()[0];
        # pacbio_header = '%s/%d/0_%d RQ=0.850' % (trimmed_header, current_read, len(read[1]));
        pacbio_header = 'S1/%d/0_%d RQ=0.850' % (current_read, len(read[1]));
        header_conversion_hash[pacbio_header] = header;
        read[0] = '%s%s' % (read[0][0], pacbio_header); ### Keep the first char of the header line.
        read[1] = re.sub("(.{500})", "\\1\n", read[1], 0, re.DOTALL);   ### Wrap the sequence line, because DALIGNER has a 9998bp line len limit.
        if (len(read) == 4):
            read[3] = re.sub("(.{500})", "\\1\n", read[3], 0, re.DOTALL);   ### Wrap the qual line, because DALIGNER has a 9998bp line len limit.
        fp_out.write('\n'.join(read) + '\n');

    sys.stderr.write('\n');
    fp_in.close();
    fp_out.close();

    return header_conversion_hash;


def add_dummy_qv(input_fastq_path, dummy_qv_val, out_fastq_path, fp_out):
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
        quals = chr(dummy_qv_val + 33) * len(read[1]);
        fp_out.write('@%s\n%s\n+\n%s\n' % (read[0][1:], read[1], quals));
        num_reads += 1;

    sys.stderr.write('\n');
    fp_in.close();

if __name__ == "__main__":
    if (len(sys.argv) < 2):
        sys.stderr.write('Various filtering methods for FASTA/FASTQ files.\n');
        sys.stderr.write('Usage:\n');
        sys.stderr.write('\theader\n');
        sys.stderr.write('\tduplicateid\n');
        sys.stderr.write('\tjoin\n');
        sys.stderr.write('\tminlen\n');
        sys.stderr.write('\tmaxlen\n');
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
        sys.stderr.write('\tqvfilter\n');
        sys.stderr.write('\tinfo\n');
        sys.stderr.write('\tcount1d2d\n');
        sys.stderr.write('\tsubsample\n');
        sys.stderr.write('\tfastq2fasta\n');
        sys.stderr.write('\tfastq2fasta2\n');
        sys.stderr.write('\tgetPairedHeaders\n');
        sys.stderr.write('\tchecknanoporepaths\n');
        sys.stderr.write('\tlength_distribution\n')
        sys.stderr.write('\t1d\n');
        sys.stderr.write('\t2d\n');
        sys.stderr.write('\treadid\n');
        sys.stderr.write('\tsubseqs\n');
        sys.stderr.write('\tmsa2fasta\n');
        sys.stderr.write('\tseparate\n');
        sys.stderr.write('\t2pacbio\n');
        sys.stderr.write('\tadddummyqv\n');

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

    elif (sys.argv[1] == 'minlen'):
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

        filter_seqs_by_length(input_fastq_path, min_length, False, out_fastq_path, fp_out);

        if (fp_out != sys.stdout):
            fp_out.close();
        exit(0);

    elif (sys.argv[1] == 'maxlen'):
        if (len(sys.argv) < 4 or len(sys.argv) > 5):
            sys.stderr.write('Extracts only sequences which have length <= max_len.\n');
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s max_len <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
            sys.stderr.write('\n');
            exit(0);

        max_length = int(sys.argv[2]);
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

        filter_seqs_by_length(input_fastq_path, max_length, True, out_fastq_path, fp_out);

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

    elif (sys.argv[1] == 'qvfilter'):
        if (len(sys.argv) < 5 or len(sys.argv) > 6):
            sys.stderr.write('Output only reads which have average base qualities above or below certain threshold.\n');
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s lte_gte threshold <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
            sys.stderr.write('\n');
            sys.stderr.write('\tlte_gte - Select which reads to output. Can be either "lt", "lte", "gt", "gte" or "eq".\n');
            sys.stderr.write('\n');
            exit(0);

        lte_gte = sys.argv[2];
        qv_threshold = int(sys.argv[3]);
        input_fastq_path = sys.argv[4];

        if ((lte_gte in ['lt', 'lte', 'gt', 'gte', 'eq']) == False):
            sys.stderr.write('ERROR: Incorrect value of the lte_gte parameter. Should be either "<" or ">".');
            exit(1);

        out_fastq_path = '';
        fp_out = sys.stdout;
        if (len(sys.argv) == 6):
            out_fastq_path = sys.argv[5];
            if (input_fastq_path == out_fastq_path):
                sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
                exit(0);
            try:
                fp_out = open(out_fastq_path, 'w');
            except Exception, e:
                sys.stderr.write(str(e));
                exit(0);

        base_quality_filter(input_fastq_path, lte_gte, qv_threshold, out_fastq_path, fp_out)

        if (fp_out != sys.stdout):
            fp_out.close();

        exit(0);

    elif (sys.argv[1] == 'info'):
        if (len(sys.argv) < 3 or len(sys.argv) > 4):
            sys.stderr.write('Tool for obtaining basic stats from FASTA/FASTQ files, such as number of sequences, total sequence length, average sequence length, etc.\n');
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s <input_fastq_file> [<input_reference_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
            sys.stderr.write('\n');
            sys.stderr.write('\t<input_reference_file> - If a reference is given, coverage is caluclated using the Lander-Waterman equation.\n');
            sys.stderr.write('\n');
            exit(0);

        input_fastq_path = sys.argv[2];
        if (len(sys.argv) == 4):
            reference_path = sys.argv[3];
            [ref_ret_string, ref_num_seqs, ref_total_seq_len, ref_average_seq_len, max_seq_len] = fastqparser.count_seq_length(reference_path);
#            sys.stdout.write('Reference info:\n');
            sys.stdout.write('(reference) Info for "%s".\n' % (reference_path));
            sys.stdout.write(ref_ret_string);
            sys.stdout.write('\n');

        [ret_string, num_seqs, total_seq_len, average_seq_len, max_seq_len] = fastqparser.count_seq_length(input_fastq_path);
#        sys.stdout.write('FASTQ info:\n');
        if (len(sys.argv) == 4):
            sys.stdout.write('(reads) Info for "%s".\n' % (input_fastq_path));
        else:
            sys.stdout.write('Info for "%s".\n' % (input_fastq_path));
        sys.stdout.write(ret_string);

        if (len(sys.argv) == 4):
            sys.stdout.write('\n');
            sys.stdout.write('Coverage: %.2fx\n' % (float(total_seq_len) / float(ref_total_seq_len)));
            sys.stdout.write('\n');

        exit(0);

    elif (sys.argv[1] == 'count1d2d'):
        if (len(sys.argv) < 3 or len(sys.argv) > 3):
            sys.stderr.write('This is not an actual filter, but counts the number of 1d or 2d reads.\n');
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s <input_fastq_file>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
            sys.stderr.write('\n');
            exit(0);

        input_fastq_path = sys.argv[2];

        count_1d2d(input_fastq_path);

        exit(0);

    elif (sys.argv[1] == 'subsample'):
        if (len(sys.argv) < 5 or len(sys.argv) > 5):
            sys.stderr.write('Subsamples a given fasta/fastq file to a given coverage. Prints the result to stdout.\n')
            sys.stderr.write('If the initail coverage is too low, prints out the original file.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <input_fasta(q)_file> <desired_coverage> <ref_genome_size>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]))
            sys.stderr.write('\n')
            exit(0);

        input_fastq_path = sys.argv[2]
        desired_coverage = int(sys.argv[3])
        ref_genome_size = int(sys.argv[4])

        subsample(input_fastq_path, desired_coverage, ref_genome_size)

        exit(0);

    elif (sys.argv[1] == 'getPairedHeaders'):
        if (len(sys.argv) != 4 ):
            sys.stderr.write('Outputs all fasta/fastq sequences for an input file, whose headers are present in a given target fasta/fastq file.\n')
            sys.stderr.write('Prints results to stdout.\n')
            sys.stderr.write('Useful when subsampling paired end reads in separate files.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <input_fasta(q)_file> <target_fasta(q)_file>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]))
            sys.stderr.write('\n')
            exit(0);

        input_fastq_path = sys.argv[2]
        target_fastq_path = sys.argv[3]

        getPaired(input_fastq_path, target_fastq_path)

        exit(0);

    elif (sys.argv[1] == 'fastq2fasta'):
        if (len(sys.argv) < 3 or len(sys.argv) > 3):
            sys.stderr.write('Outputs a given fastq file in fasta format. Replaces the ":" characters with blank spaces.\n')
            sys.stderr.write('Additionally, replaces colons in header with spaces (relevant for running L&S pipeline).\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <input_fastq_file>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]))
            sys.stderr.write('\n')
            exit(0);

        input_fastq_path = sys.argv[2]

        fastq2fasta(input_fastq_path, True)

        exit(0);

    elif (sys.argv[1] == 'fastq2fasta2'):
        if (len(sys.argv) < 3 or len(sys.argv) > 3):
            sys.stderr.write('Outputs a given fastq file in fasta format. Does not replace the ":" characters with blank spaces (unlike fastq2fasta).\n')
            sys.stderr.write('Additionally, replaces colons in header with spaces (relevant for running L&S pipeline).\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <input_fastq_file>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]))
            sys.stderr.write('\n')
            exit(0);

        input_fastq_path = sys.argv[2]

        fastq2fasta(input_fastq_path, False)

        exit(0);

    elif (sys.argv[1] == 'checknanoporepaths'):
        if (len(sys.argv) < 4 or len(sys.argv) > 4):
            sys.stderr.write('This is not an actual filter, but removes extra spaces from nanopore paths.\n')
            sys.stderr.write('So that reads file can be used with nanopolish.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <input_fastq_file> <fast5_root_path>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]))
            sys.stderr.write('\n')
            exit(0)

        input_fastq_path = sys.argv[2]
        fast5_root_path = sys.argv[3]

        check_nanopore_paths(input_fastq_path, fast5_root_path)

        exit(0)

    elif (sys.argv[1] == 'length_distribution'):
        if (len(sys.argv) < 3 or len(sys.argv) > 3):
            sys.stderr.write('Outputs the distribution of reads in a fasta/fastq file.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <input_file>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]))
            sys.stderr.write('\n')
            exit(0)

        input_path = sys.argv[2]
        length_distribution(input_path)

        exit(0);



    elif (sys.argv[1] == '1d' or sys.argv[1] == '2d'):
        if (len(sys.argv) < 3 or len(sys.argv) > 4):
            sys.stderr.write('Extracts only reads with typical %s (nanopore) headers (those containing either "1d", "template" or "complement"; "2d" or "twodir" in their header).\n' % (sys.argv[1]));
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
            sys.stderr.write('\n');
            exit(0);

        # header_patterns_path = sys.argv[2];
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

        if (sys.argv[1] == '1d'):
            filter_seqs_by_header_list(input_fastq_path, ['1d', 'template', 'complement'], out_fastq_path, fp_out);
        else:
            filter_seqs_by_header_list(input_fastq_path, ['2d', 'twodir'], out_fastq_path, fp_out);

        if (fp_out != sys.stdout):
            fp_out.close();
        exit(0);

    elif (sys.argv[1] == 'readid'):
        if (len(sys.argv) < 4 or len(sys.argv) > 5):
            sys.stderr.write('Takes a FASTA/FASTQ file and a file containing seq IDs (0-offset). Extracts all seqs with given IDs.\n');
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s <input_id_file> <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
            sys.stderr.write('\n');
            exit(0);

        read_id_path = sys.argv[2];
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

        filter_seqs_by_read_id(input_fastq_path, read_id_path, out_fastq_path, fp_out);

        if (fp_out != sys.stdout):
            fp_out.close();
        exit(0);

    elif (sys.argv[1] == 'subseqs'):
        if (len(sys.argv) < 5 or len(sys.argv) > 6):
            sys.stderr.write('Extracts bases from all sequences in a FASTA file between specified coordinates. End coordinate is not inclusive.\n');
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s <input_fastq_file> start end [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
            sys.stderr.write('\n');
            exit(0);

        input_fastq_path = sys.argv[2];
        start_coord = int(sys.argv[3]);
        end_coord = int(sys.argv[4]);

        out_fastq_path = '';
        fp_out = sys.stdout;
        if (len(sys.argv) == 6):
            out_fastq_path = sys.argv[5];
            if (input_fastq_path == out_fastq_path):
                sys.stderr.write('ERROR: Output and input files are the same! Exiting.\n');
                exit(0);
            try:
                fp_out = open(out_fastq_path, 'w');
            except Exception, e:
                sys.stderr.write(str(e));
                exit(0);

        extract_subseqs(input_fastq_path, start_coord, end_coord, out_fastq_path, fp_out);

        if (fp_out != sys.stdout):
            fp_out.close();
        exit(0);

    elif (sys.argv[1] == 'msa2fasta'):
        if (len(sys.argv) < 3 or len(sys.argv) > 4):
            sys.stderr.write('Takes a multifasta file of multiple sequence alignments, and produces a majority vote consensus. Only [ACTG] bases are considered for output. Gaps in the input file are denoted by either a \'.\' or a \'-\'.\n');
            sys.stderr.write('Usage:\n')
            sys.stderr.write('\t%s %s <input_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]))
            sys.stderr.write('\n')
            exit(0)

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

        msa2fasta(input_fastq_path, fp_out)

        if (fp_out != sys.stdout):
            fp_out.close();
        exit(0);        

    elif (sys.argv[1] == 'separate'):
        if (len(sys.argv) < 4 or len(sys.argv) > 4 ):
            sys.stderr.write('Separate a FASTA/FASTQ file into individual sequence files. File names will be numbers assigned to sequences sequentially.\n');
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s <input_fastq_file> <out_folder>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
            sys.stderr.write('\n');
            exit(0);

        input_fastq_path = sys.argv[2];
        out_folder = sys.argv[3];

        separate_seqs(input_fastq_path, out_folder);

    elif (sys.argv[1] == '2pacbio'):
        if (len(sys.argv) < 4 or len(sys.argv) > 4):
            sys.stderr.write('Makes all headers in a FASTQ file unique by adding a number at the end of the header. Also, splits headers up to the first whitespace.\n');
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s <input_fastq_file> <out_converted_fastq_file>\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
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
        fp_out.close();

        convert_reads_to_pacbio_format(input_fastq_path, out_fastq_path)

        # if (fp_out != sys.stdout):
        #     fp_out.close();

        exit(0);

    elif (sys.argv[1] == 'adddummyqv'):
        if (len(sys.argv) < 4 or len(sys.argv) > 5):
            sys.stderr.write('Takes an input FASTA/FASTQ file and converts it to a FASTQ format if necessary, while replacing all QVs with a given qv_val value.\n');
            sys.stderr.write('Usage:\n');
            sys.stderr.write('\t%s %s qv_val <input_fastq_file> [<out_filtered_fastq_file>]\n' % (os.path.basename(sys.argv[0]), sys.argv[1]));
            sys.stderr.write('\n');
            exit(0);

        dummy_qv_val = int(sys.argv[2]);
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

        add_dummy_qv(input_fastq_path, dummy_qv_val, out_fastq_path, fp_out);

        if (fp_out != sys.stdout):
            fp_out.close();
        exit(0);


    else:
        sys.stderr.write('ERROR: Unknown subcommand!\n');
        exit(0);
