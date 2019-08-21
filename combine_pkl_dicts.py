'''
Takes as input a list of pickle files containing dictionaries,
imports them, combines them into a single dict,
then writes a new pickle

Note that there should be no common key
across the various provided dicts.  If there are,
some entries will overwrite others.
'''

import cPickle as pkl
import argparse
import sys
import subprocess



def parse_arguments(input_args):

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--pkl_paths",
		type = str,
		help = ("Comma separated paths " +
				"pickled dict files"),
		required = True)

	parser.add_argument(
		"--outdir",
		type = str,
		help = "Path to output directory",
		required = True)

	parser.add_argument(
		"--prefix",
		type = str,
		help = "Output file prefix",
		required = True)

	args = parser.parse_args(input_args)

	return(args)


def main(args):

	subprocess.call("mkdir -p ", args.outdir, shell = True)

	pkl_path_list = args.cds_ins_pkl.split(",")

	print "Reading pickle files"

	dict_list = [pkl.load(open(i, "rb")) for i in pkl_path_list]

	print "Combining list of dicts into single dict"

	combined_dict = dict(chain(*map(dict.items, dict_list)))

	print "Writing output pickle file"

	pkl.dump(combined_dict, 
		open(args.outdir + "/" + args.prefix + "_combined_dict.pkl", "wb"))

	print "Finished"



if __name__ == "__main__":

	output_args <- parse_arguments(sys.argv[1:])
	main(output_args)