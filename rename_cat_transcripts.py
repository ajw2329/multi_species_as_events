'''
Rename CAT transcripts to exactly 
match original IDs when possible
'''


import argparse
from splice_lib import splice_lib
import sys
from collections import Counter
import copy



def check_ids(
		ref_transcript_dict,
		cat_transcript_dict,
		sep = "-"):

	renamed_appearance = [] ## count # of times renamed id appears - 
							## in case of multiple tx versions in CAT

	renamed_cat_transcript_dict = {}

	for transcript, transcript_val in cat_transcript_dict.iteritems():

		for substr in transcript.split("-"):

			if substr in ref_transcript_dict:

				renamed_appearance.append(substr)

				renamed_cat_transcript_dict[substr] = copy.deepcopy(transcript_val)

	count_renamed_appearance = Counter(renamed_appearance)

	print len(count_renamed_appearance)

	non_unique_count = 0

	for tx, count in count_renamed_appearance.iteritems():

		if count > 1:

			non_unique_count += 1

	print "There are", non_unique_count, "non-unique ids."

	return renamed_cat_transcript_dict


def parse_arguments(input_args):

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--reference_gtf",
		type = str,
		help = "Path to reference tx gtf",
		required = True)

	parser.add_argument(
		"--cat_gtf",
		type = str,
		help = "Path to CAT-produced non-ref tx gtf",
		required = True)

	parser.add_argument(
		"--outdir",
		type = str,
		help = "Path to outdir",
		required = True)

	parser.add_argument(
		"--prefix",
		type = str,
		help = "Prefix for output file name",
		required = True)


	args = parser.parse_args()

	return args


def main(args):

	ref_transcript_dict = splice_lib.generate_standard_transcript_dict(
		args.reference_gtf)

	cat_transcript_dict = splice_lib.generate_standard_transcript_dict(
		args.cat_gtf)

	renamed_cat_dict = check_ids(
		ref_transcript_dict,
		cat_transcript_dict,
		sep = "-")

	splice_lib.output_transcript_gtf(
		renamed_cat_dict,
		args.outdir,
		name = args.prefix + "_renamed_to_ref")


if __name__ == "__main__":

	args = parse_arguments(sys.argv[1:])
	main(args)

