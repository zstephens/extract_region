import argparse
import pysam
import re
import os

MATCH_CHAR = 'MX='
REF_CHAR   = 'MX=D'
READ_CHAR  = 'MX=IS'

def get_ref_2_read_dict(read_dat):
	cigar   = read_dat[3]
	letters = re.split(r"\d+",cigar)[1:]
	numbers = [int(n) for n in re.findall(r"\d+",cigar)]
	[adj, radj] = [0,0]
	refpos_to_readpos = {}
	for i in range(len(letters)):
		if letters[i] in MATCH_CHAR:
			for j in range(numbers[i]):
				refpos_to_readpos[read_dat[2]+adj+j] = radj+j
		if letters[i] in REF_CHAR:
			for j in range(numbers[i]):
				refpos_to_readpos[read_dat[2]+adj+j] = radj+j
			adj += numbers[i]
		if letters[i] in READ_CHAR:
			radj += numbers[i]
	return refpos_to_readpos

def exists_and_is_nonzero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

def main(raw_args=None):
	parser = argparse.ArgumentParser(description='extract_region.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',  type=str, required=True,                 metavar='input.bam',     help="* Input BAM")
	parser.add_argument('-o',  type=str, required=True,                 metavar='output.fa',     help="* Output sequences")
	parser.add_argument('-c',  type=str, required=True,                 metavar='chr:start-end', help="* Reference coordinates")
	parser.add_argument('-b',  type=int, required=False, default=50000, metavar='50000',         help="Buffer size for extracting reads surrounding specified region")
	parser.add_argument('-m',  type=int, required=False, default=3,     metavar='3',             help="Minimum read MAPQ")
	parser.add_argument('--keep-suppl',  required=False, default=False, action='store_true',     help='Do not discard supplementary alignments')
	args = parser.parse_args()

	INPUT_ALN    = args.i
	OUT_FASTA    = args.o
	IN_REGION    = args.c
	SAMVIEW_BUFF = args.b
	MAPQ_THRESH  = args.m
	KEEP_SUPPL   = args.keep_suppl
	#
	splt  = IN_REGION.split(':')
	splt2 = splt[1].split('-')
	(MY_CHR, MY_START, MY_END) = (splt[0], int(splt2[0]), int(splt2[1])+1)	# start-end are 1-indexed and inclusive
	samstart  = max([0, MY_START - SAMVIEW_BUFF])
	samend    = MY_END + SAMVIEW_BUFF
	samregion = MY_CHR + ':' + str(samstart) + '-' + str(samend)
	#
	if exists_and_is_nonzero(INPUT_ALN) == False:
		print('Error reading -i input file')
		exit(1)
	if INPUT_ALN[-4:].lower() == '.sam':
		samfile = pysam.AlignmentFile(INPUT_ALN, "r")
	elif INPUT_ALN[-4:].lower() == '.bam':
		samfile = pysam.AlignmentFile(INPUT_ALN, "rb")
	elif INPUT_ALN[-4:].lower() == '.cram':
		samfile = pysam.AlignmentFile(INPUT_ALN, "rc")
	else:
		print('Error: unknown -i file type')
		exit(1)
	#
	read_list = []
	try:
		for aln in samfile.fetch(region=samregion):
			splt = str(aln).split('\t')
			if KEEP_SUPPL == False and int(splt[1])&2048:	# discard supplementary alignments
				continue
			if int(splt[4]) < MAPQ_THRESH:					# discard low mapq reads
				continue
			orientation = 'FWD'
			if int(splt[1])&16:
				orientation = 'REV'
			# rname, chr, pos, cigar, seq, qual, orientation
			read_list.append([splt[0], MY_CHR, int(splt[3]), splt[5], splt[9], splt[10], orientation])	
	except ValueError:
		print('Error:', MY_CHR, 'was not found in input BAM')
		samfile.close()
		exit(1)
	samfile.close()
	#
	seq_lens  = []
	reads_skipped = {'s_only':0, 'e_only':0, 'neither':0}
	#
	fa_out = open(OUT_FASTA, 'w')
	for i in range(len(read_list)):
		refpos_to_readpos = get_ref_2_read_dict(read_list[i])
		if MY_START in refpos_to_readpos and MY_END in refpos_to_readpos:
			read_sequence = read_list[i][4][refpos_to_readpos[MY_START]:refpos_to_readpos[MY_END]]
			seq_lens.append(len(read_sequence))
			fa_out.write('>' + read_list[i][0] + '\n' + read_sequence + '\n')
		elif MY_START in refpos_to_readpos and MY_END not in refpos_to_readpos:
			reads_skipped['s_only'] += 1
		elif MY_START not in refpos_to_readpos and MY_END in refpos_to_readpos:
			reads_skipped['e_only'] += 1
		else:
			reads_skipped['neither'] += 1
	fa_out.close()
	#
	print('sequences written:', len(seq_lens))
	print('reads skipped:')
	print(' - spanned start coord but not end:', reads_skipped['s_only'])
	print(' - spanned end coord but not start:', reads_skipped['e_only'])
	print(' - spanned neither start nor end:  ', reads_skipped['neither'])

if __name__ == '__main__':
	main()
