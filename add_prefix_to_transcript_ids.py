'''
Reads in a transcript GTF, and a prefix

Writes a new transcript GTF with transcript ids
containing the prefix
'''

import argparse
from splice_lib import splice_lib
import sys
import subprocess




def parse_arguments(input_args):

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--transcript_gtf",
		type = str,
		help = "Path to transcript GTF",
		required = True)

	parser.add_argument(
		"--prefix",
		type = str,
		help = "Prefix for transcript ids and outfile",
		required = True)

	parser.add_argument(
		"--outdir",
		type = str,
		help = "Path to output directory",
		required = True)

	parser.add_argument(
		"--outfile_suffix",
		type = str,
		help = "Optional outfile suffix to override 'splice_lib_transcripts'",
		default = "splice_lib_transcripts"
		)

	args = parser.parse_args(input_args)

	return(args)


def main(args):

	subprocess.call("mkdir -p " + args.outdir, shell = True)

	print "Importing standard transcript dict"

	transcript_dict = splice_lib.generate_standard_transcript_dict(args.transcript_gtf)

	print "Writing prefixed standard_transcript_dict to GTF"
	splice_lib.output_transcript_gtf(
        transcript_dict, 
        args.outdir, 
        name = args.prefix + "_" + args.outfile_suffix,
        transcript_prefix = args.prefix)


if __name__ == "__main__":

	output_args = parse_arguments(sys.argv[1:])
	main(output_args)