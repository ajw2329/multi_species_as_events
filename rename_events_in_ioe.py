'''
Rename event ids in splice_lib_events.ioe files

Intended for situations where new consensus event
ids are created for orthologous events in 
multiple species.  

Takes as input an ioe file and final_original_name.tsv
(from multi_species_events_standalone.py) and writes
new ioe file with the consensus event id from
final_original_name.tsv replacing the original
'''

import argparse
import sys
import subprocess


def create_name_dict(
		final_original_name,
		original_index):

	name_dict = {}

	with open(final_original_name, 'r') as file:

		next(file) # skip header

		for line in file:

			entry = line.split()

			old_name = entry[original_index]
			new_name = entry[0]

			name_dict[old_name] = new_name

	return name_dict

	

def add_transcripts_to_event_dict(
        ioe_file, 
        name_dict,
		outdir,
		prefix,
		remove_unmatched = False):

	'''
        Creates pseudo event and transcript dicts 
        to run event_nmd_nsd_status and output_table 
        functions without having passed the full 
        dicts as input to main.  Useful if 
        running script as standalone.
	'''

	outfile = open(outdir + "/" + prefix + "_splice_lib_events.ioe", 'w')

	with open(ioe_file, 'r') as file:

		header = next(file) # skip header
		outfile.write(header)

		for line in file:

			entry = line.strip().split()
			event = entry[2].split(";")[1]

			included_form_transcripts = set(entry[3].split(","))
			all_transcripts = entry[4].split(",")
			excluded_form_transcripts = set(all_transcripts) - included_form_transcripts

			new_name = name_dict.get(event)

			if new_name:

				outname = entry[2].split(";")[0] + ";" + new_name

			else:

				if not remove_unmatched:

					outname = entry[2]

				else:

					continue

			out_entry = "\t".join([
				entry[0],
				entry[1],
				outname,
				entry[3],
				entry[4]])

			outfile.write(out_entry + "\n")

	outfile.close()



def parse_arguments(input_args):

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--ioe_file",
		type = str,
		help = "Path to input ioe file",
		required = True)

	parser.add_argument(
		"--final_original_name",
		type = str,
		help = ("Path to final_original_name.tsv " +
			    "which contains new event_ids"),
		required = True)

	parser.add_argument(
		"--original_index",
		type = int,
		help = "Index of original event id in final_original_name",
		required = True)

	parser.add_argument(
		"--outdir",
		type = str,
		help = "Path to output directory",
		required = True)

	parser.add_argument(
		"--prefix",
		type = str,
		help = "Prefix for output file",
		default = "renamed")

	parser.add_argument(
		"--remove_unmatched",
		action = "store_true",
		help = "Remove events that don't have a match in final_original_name")


	output_args = parser.parse_args(input_args)

	return output_args


def main(args):

	name_dict = create_name_dict(
					args.final_original_name,
					args.original_index)

	add_transcripts_to_event_dict(
        args.ioe_file, 
        name_dict,
		args.outdir,
		args.prefix,
		args.remove_unmatched)

if __name__ == "__main__":

	output_args = parse_arguments(sys.argv[1:])

	main(output_args)

