'''
Intended to parse spladder events lifted from GRCh38 to Rhemac8 or vice versa

Takes as input a gtf file and returns a two-field tsv with the matched ID of each species in each field
'''

import re
import sys
from splice_lib import splice_lib
import argparse
import copy
import cPickle as pkl
import subprocess


def index_genes_by_junctions(gene_dict):
	'''
		Returns dictionary where splice junctions (specified by chrom, strand, and inner coordinates of flanking exons) are keys which point to lists of matching genes

		Takes the transcript dictionary as input.

		This facilitates the matching of CDS without needing prior gene information.

	'''

	junction_dict = {}

	for transcript in gene_dict:

		gene_dict[transcript]["all_jn"] = splice_lib.get_junctions(gene_dict[transcript]["exons"])

		for junction in gene_dict[transcript]["all_jn"]:

			junction_key = gene_dict[transcript]["chrom"] + "_" + "_".join(map(str,junction)) + "_" + gene_dict[transcript]["strand"]

			junction_dict.setdefault(junction_key, set()).add(gene_dict[transcript]["gene"])

	return junction_dict
	

def annotate_custom_transcripts(
		standard_transcript_dict, 
		gene_junction_dict, 
		transcripts = None, 
		prev_length = None):

	if transcripts is None:

		transcripts = standard_transcript_dict.keys()
		prev_length = len(transcripts)

	expanded_gene_junction_dict = copy.deepcopy(gene_junction_dict)

	unmatched_tx = set()
	matched_tx = set()

	for transcript in transcripts:

		transcript_val = standard_transcript_dict[transcript]

		transcript_val["gene_ambiguous_in_annotation"] = False

		transcript_junctions = splice_lib.get_junctions(transcript_val["exons"])

		junction_keys = [
			transcript_val["chrom"] + 
			"_" + 
			"_".join(map(str, junction)) + 
			"_" + 
			transcript_val["strand"] for 
			junction in transcript_junctions]

		compatible_genes = []

		for junction_key in junction_keys:

			if junction_key in gene_junction_dict:

				if len(gene_junction_dict[junction_key]) > 1:

					transcript_val["gene_ambiguous_in_annotation"] = True

				compatible_genes.extend(gene_junction_dict[junction_key])

		compatible_genes = list(set(compatible_genes))

		if len(compatible_genes) > 0:

			matched_tx.add(transcript)

			for junction_key in junction_keys:

				if junction_key in gene_junction_dict:

					expanded_gene_junction_dict[junction_key] = list(
						set(compatible_genes + 
							list(gene_junction_dict[junction_key])))

				else:

					expanded_gene_junction_dict[junction_key] = compatible_genes

			transcript_val["gene"] = ",".join(compatible_genes)


		else:

			unmatched_tx.add(transcript) 

	len_unmatched = len(unmatched_tx)

	if len_unmatched < prev_length:

		return annotate_custom_transcripts(
					standard_transcript_dict, 
					expanded_gene_junction_dict, 
					transcripts = unmatched_tx, 
					prev_length = len_unmatched)

	else:

		return expanded_gene_junction_dict, unmatched_tx






def find_matches(
		event_dict, 
		event_keys, 
		gene_junction_dict):

	matches = {}

	unmatched = []


	for event in event_keys:

		matched = False

		junction_list = [
			list(i) for i in 
			set(tuple(i) for i in 
				(event_dict[event]["included_junctions"] + 
				 event_dict[event]["excluded_junctions"]))]

		for junction in junction_list:

			junction_key = (
				event_dict[event]["chrom"] + 
				"_" + 
				"_".join(map(str,junction)) + 
				"_" + 
				event_dict[event]["strand"])


			if junction_key in gene_junction_dict:

				matched = True

				for match in gene_junction_dict[junction_key]:

					if match not in matches:

						matches[match] = []

					if event not in matches[match]:
						matches[match].append(event)

		if matched == False:

			unmatched.append(event)

	return [matches, unmatched]



def generate_transcript_chrom_dict(chrom_transcript):

	transcript_chrom_dict = {}

	with open(chrom_transcript) as file:

		for line in file:

			entry = line.split()
			chrom = entry[0]
			transcript_id = re.sub('[";]', "", entry[1]).split(".")[0]

			transcript_chrom_dict[transcript_id] = chrom


	return transcript_chrom_dict


def replace_chrom(genes_dict, transcript_chrom_dict):

	for transcript in genes_dict:

		if transcript.split(".")[0] in transcript_chrom_dict:

			genes_dict[transcript]["chrom"] = transcript_chrom_dict[transcript.split(".")[0]]

		else:
			print "Missing", transcript


def generate_transcript_gene_dict(transcript_gene):

	transcript_gene_dict = {}

	with open(transcript_gene) as file:

		for line in file:

			entry = line.split()
			transcript_id = re.sub('[";]', "", entry[0]).split(".")[0].split(".")[0]
			gene_id = re.sub('[";]', "", entry[1]).split(".")[0]

			transcript_gene_dict[transcript_id] = gene_id

	return transcript_gene_dict


def add_gene_ids(genes_dict, transcript_gene_dict):

	for transcript in genes_dict:

		if transcript.split(".")[0] in transcript_gene_dict:

			genes_dict[transcript]["gene_id"] = transcript_gene_dict[transcript.split(".")[0]]

		else:

			print "Missing", transcript


def main(args, event_dict = None):

	parser = argparse.ArgumentParser()

	parser.add_argument(
		'--event_gtf', 
		type = str, 
		help = "GTF file of species 1 reference events")

	parser.add_argument(
		'--gene_gtf', 
		type = str, 
		help = "GTF file of full genes - GENCODE format expected", 
		required = True)

	parser.add_argument(
		'--transcript_gtf', 
		type = str, 
		help = ("Optional GTF file of transcripts e.g. " + 
			    "generated by StringTie - allows more " + 
			    "comprehensive event assignment as many-junction " + 
			    "transcripts are more likely to have " + 
			    "overlap with annotated genes"))

	parser.add_argument(
		'--prefix', 
		type = str, 
		help = "Prefix to render output file unique", 
		required = True)

	parser.add_argument(
		'--outdir', 
		type = str, 
		help = "Path to output directory")

	parser.add_argument(
		'--chrom_transcript', 
		type = str, 
		help = ("File associating chromosome names " + 
			    "and transcript ids - useful for " + 
			    "changing chrom IDs (e.g. from ensembl to ucsc"))

	parser.add_argument(
		'--transcript_gene', 
		type = str, 
		help = ("File associating genes to transcript ids " + 
			    "- useful when gene IDs (symbols) not " + 
			    "contained in original gtf"))

	parser.add_argument(
		'--suppress_output', 
		action = "store_true", 
		help = "Suppress output files")

	args = parser.parse_args(args)

	event_gtf = args.event_gtf
	gene_gtf = args.gene_gtf
	transcript_gtf = args.transcript_gtf
	outdir = args.outdir
	prefix = args.prefix
	chrom_transcript = args.chrom_transcript
	transcript_gene = args.transcript_gene
	suppress_output = args.suppress_output


	if event_gtf is not None and event_dict is not None:

		sys.exit("Can't provide both event_dict and event_gtf - assign_events_to_genes does not know which to use! Exiting . . . ")

	if not suppress_output and not outdir:

		sys.exit("--outdir must be specified if --suppress_output is not. Exiting . . . ")

	
	if outdir:

		subprocess.call("mkdir -p " + outdir, shell = True)


	genes_dict = splice_lib.generate_standard_transcript_dict(gene_gtf)
	splice_lib.sort_transcript_dict_exons(genes_dict)

	if chrom_transcript != None:

		transcript_chrom_dict = generate_transcript_chrom_dict(chrom_transcript)
		replace_chrom(genes_dict, transcript_chrom_dict)

	if transcript_gene != None:

		transcript_gene_dict = generate_transcript_gene_dict(transcript_gene)
		add_gene_ids(genes_dict, transcript_gene_dict)

	if not suppress_output:

		with open(outdir + '/' + prefix + '_genes_annotated_transcripts.tsv', 'w') as file:

			for transcript in genes_dict:

				file.write(genes_dict[transcript]["gene"] + "\t" + transcript + "\n")


	gene_junction_dict = index_genes_by_junctions(genes_dict)

	print "gene_junction_dict_len", len(gene_junction_dict)

	multiply_assigned_junction_dict = {}

	multiply_assigned_junction_counter = 0

	for junction in gene_junction_dict:

		if len(gene_junction_dict[junction]) > 1:

			multiply_assigned_junction_counter += 1

			multiply_assigned_junction_dict[junction] = gene_junction_dict[junction]

	print "multiply assigned junction count", multiply_assigned_junction_counter

	if not suppress_output:

		pkl.dump(multiply_assigned_junction_dict, open(outdir + "/multiply_assigned_junctions.pkl", "wb"))

	if transcript_gtf != None:

		standard_transcript_dict = splice_lib.generate_standard_transcript_dict(transcript_gtf)
		splice_lib.sort_transcript_dict_exons(standard_transcript_dict)
		expanded_gene_junction_dict, unmatched_tx = annotate_custom_transcripts(standard_transcript_dict, gene_junction_dict)

		if not suppress_output:

			with open(outdir + '/' + prefix + '_genes_custom_transcripts.tsv', 'w') as file:

				for transcript in standard_transcript_dict:

					file.write(standard_transcript_dict[transcript]["gene"] + "\t" + transcript + "\n")

	else:

		standard_transcript_dict = None

	if event_gtf is not None or event_dict is not None:

		if event_gtf is not None:

			events_dict = splice_lib.generate_standard_event_dict(event_gtf)
			splice_lib.complete_event_dict(events_dict)

		elif event_dict is not None:

			events_dict = event_dict

		matched_results = find_matches(events_dict, events_dict.keys(), gene_junction_dict)

		matches_dict = matched_results[0]

		unmatched = matched_results[1]

		### now create new gene junction dict based on transcript_gtf

		if transcript_gtf != None:

			supplemented_gene_junction_dict = index_genes_by_junctions(standard_transcript_dict)

			new_matched_results = find_matches(events_dict, unmatched, expanded_gene_junction_dict)

			new_match_dict = new_matched_results[0]

			unmatched = new_matched_results[1]

			for gene in new_match_dict:

				if gene in matches_dict:

					matches_dict[gene].extend(new_match_dict[gene])

				else:

					matches_dict[gene] = copy.deepcopy(new_match_dict[gene])

		if not suppress_output:

			with open(outdir + '/' + prefix + '_genes_events.tsv', 'w') as file:

				for gene in matches_dict:

					for event in matches_dict[gene]:

						file.write(gene + "\t" + event + "\n")

			with open(outdir + '/' + prefix + '_orphan_events.tsv', 'w') as file:

				for event in unmatched:

					file.write(event + "\n")

	else:
		matches_dict = None


	return matches_dict, genes_dict, standard_transcript_dict


if __name__ == '__main__':

	main(sys.argv[1:])
