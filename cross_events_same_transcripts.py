'''
Program to compare pairwise events generated from a single annotation that has been mapped to two species (the original, and a closely related species)

Primary use case as of 10/27/2017: Compare events generated from gencode v24 (GRCh38) and gencode v24 mapped to RheMac8 by Ian Fiddes (Haussler Lab) and
generate a table of orthologous events
'''

import argparse
import pickle
import sys
from splice_lib import splice_lib


##TODO make compatible with assess_junction_support

def import_transcript_matching_dict(
	transcript_table,
	transcript_dict_1):
	'''
		Imports simple two-field tsv in 
		which annotation 1 transcript IDs 
		are in first field (i.e. IDs 
		corresponding to transcripts in 
		transcript_gtf_1) and annotation 
		2 transcript IDs are in the second field.  
		Returns dict whose keys are annotation 1 
		IDs and whose values are annotation 2 IDs.
	'''

	transcript_match_dict = {}

	if transcript_table is not None:

		with open(transcript_table, 'r') as file:

			for line in file:

				entry = line.strip().split()

				transcript_1 = entry[0]
				transcript_2 = entry[1]

				if transcript_1 not in transcript_match_dict:

					transcript_match_dict[transcript_1] = transcript_2

				else:
					print ("Warning - " +
					       transcript_1 +
					       " appears twice in input file - " + 
					       "implies non-one-to-one transcript " + 
					       "relationship.  Skipping repeated instance.")

	else:

		for transcript in transcript_dict_1:

			transcript_match_dict.setdefault(
				transcript, 
				transcript)


	return transcript_match_dict




def make_transcript_event_dict(event_dict, form):
	'''
		Makes transcript-centric event dict in 
		which keys are transcript IDs and values 
		are lists of event IDs.  It is called with 
		both a splice_lib standard event dict and 
		the event form in question (i.e. either 
		"included" or "excluded") as unnamed arguments.  
		Returns said dictionary.
	'''

	transcript_centric_event_dict = {}

	for event in event_dict:

		for transcript in event_dict[event][form + "_form_transcripts"]:

			if transcript not in transcript_centric_event_dict:

				transcript_centric_event_dict[transcript] = []

			if event not in transcript_centric_event_dict[transcript]:

				transcript_centric_event_dict[transcript].append(event)

	return transcript_centric_event_dict



def find_matches(
	splice_dict_1, 
	splice_dict_2, 
	transcript_match_dict, 
	transcript_dict_1, 
	transcript_dict_2, 
	species_1_transcript_centric_included_forms, 
	species_1_transcript_centric_excluded_forms, 
	species_2_transcript_centric_included_forms, 
	species_2_transcript_centric_excluded_forms):

	'''Identifies events consisting of consistently ordered exons in at least one pair of orthologous transcripts
	'''

	event_matches = {}

	##loop through included and excluded event forms

	non_concordant_transcripts = set()

	for form in ["included", "excluded"]:

		##loop through transcripts in transcript-indexed event dict (for particular event form)

		for transcript, transcript_val in eval(
			'species_1_transcript_centric_' + 
			form + 
			'_forms').iteritems():

			##Look through all events in this particular transcript 

			for event in transcript_val:

				transcript_junction_indices_1 = []  ##Begin collecting junction indices (i.e. this is a proxy for relative position of participating exons)

				event_entry = splice_dict_1.get(event)
				event_type = event_entry["event_type"]

				for junction in event_entry[form + "_junctions"]:  ##loop through junctions in event

					if junction in transcript_dict_1[transcript]["junctions"]:  ##make sure corresponding junction exists in transcript dict

						transcript_junction_indices_1.append(
							transcript_dict_1[transcript]["junctions"].index(junction)) ##Append the junction index (relative position) in the transcript junction list

				if transcript in transcript_match_dict: ##find the corresponding transcript in the match dict

					match = transcript_match_dict[transcript]

					if (match in 
						eval(
							'species_2_transcript_centric_' + 
							form + 
							'_forms')): 

						if (len(transcript_dict_1[transcript]["junctions"]) == 
							len(transcript_dict_2[match]["junctions"])):

							for p_event in eval(
								'species_2_transcript_centric_' + 
								form + 
								'_forms[match]'):  ##repeat the process for the second annotation

								transcript_junction_indices_2 = []

								p_event_entry = splice_dict_2.get(p_event)
								p_event_type = p_event_entry["event_type"]

								for junction in p_event_entry[form + "_junctions"]:

									if junction in transcript_dict_2[match]["junctions"]:

										transcript_junction_indices_2.append(
											transcript_dict_2[match]["junctions"].index(junction))

								if ((transcript_junction_indices_1 == 
								    transcript_junction_indices_2) and 
								    (event_type == p_event_type)): ###make sure the junctions have the same relative position for the putatively matched event, and that the events have the same type 

									((event_matches.setdefault(event, {})
										).setdefault(form + "_matches", [])
											).append(p_event)

						else:

							print ("Warning: Transcript " +
							       transcript +
							       " not exon count concordant. " + 
							       "Skipping . . . ")

							non_concordant_transcripts.add(transcript)


	final_matches = {}

	inverted_matches = {}

	more_than_one_match_counter = 0

	print "There are " + str(len(event_matches)) + " events with possible matches"

	for event, event_val in event_matches.iteritems():

		if "included_matches" not in event_val:

			matches = list(set(event_val["excluded_matches"]))
			support = "excluded"

		elif "excluded_matches" not in event_val:

			matches = list(set(event_val["included_matches"]))
			support = "included"

		else:

			matches = list(set(event_val["included_matches"]) & 
			   set(event_val["excluded_matches"]))
			support = "both"

		if len(matches) == 1:

			final_matches[event] = {"match": matches[0], "support": support}

			inverted_matches.setdefault(matches[0], set()).add(event)

		elif len(matches) > 1:

			more_than_one_match_counter += 1

	for match, events in inverted_matches.iteritems():

		if len(events) > 1:

			more_than_one_match_counter += 1

			for event in events:

				del final_matches[event]


	return final_matches



def main():

	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--event_gtf_1", 
		type = str, 
		help = "GTF of first species events")

	parser.add_argument(
		"--event_gtf_2", 
		type = str, 
		help = "GTF of second species events")

	parser.add_argument(
		"--event_ioe_1",
		type = str,
		help = "IOE file of first species event-transcript associations")

	parser.add_argument(
		"--event_ioe_2",
		type = str,
		help = "IOE file of second species event-transcript associations")

	parser.add_argument(
		"--transcript_gtf_1", 
		type = str, 
		help = "GTF file for full transcript annotations species 1")

	parser.add_argument(
		"--transcript_gtf_2", 
		type = str, 
		help = "GTF file for full transcript annotations species 2")

	parser.add_argument(
		"--transcript_table", 
		type = str, 
		help = ("Two field tsv relating transcript IDs - " + 
			    "transcripts for splice_dict_1 in first field, " + 
			    "2 in second field. If not given, it will be assumed " +
			    "that transcripts have the same IDs in both species " + 
			    "already."))

	parser.add_argument(
		"--outdir", 
		type = str, 
		help = "Path to output directory", 
		required = True)

	args = parser.parse_args()



	splice_dict_1 = args.event_gtf_1
	splice_dict_2 = args.event_gtf_2
	transcript_table = args.transcript_table
	outdir = args.outdir
	transcript_gtf_1 = args.transcript_gtf_1
	transcript_gtf_2 = args.transcript_gtf_2



	splice_dict_1 = splice_lib.generate_standard_event_dict(
									args.event_gtf_1)

	print splice_dict_1.keys()[0:10]

	print "CO.0003623 splice dict 1", ("CO.0003623" in splice_dict_1)

	splice_lib.complete_event_dict(
				splice_dict_1)

	print splice_dict_1.keys()[0:10]
	print "CO.0003623 splice dict 1", ("CO.0003623" in splice_dict_1)	

	splice_lib.add_transcripts_to_event_dict(
		args.event_ioe_1, 
		splice_dict_1)

	splice_dict_2 = splice_lib.generate_standard_event_dict(
									args.event_gtf_2)

	print splice_dict_2.keys()[0:10]	

	print "CO.0003623 splice dict 2", ("CO.0003623" in splice_dict_2)	

	splice_lib.complete_event_dict(
				splice_dict_2)

	print splice_dict_2.keys()[0:10]

	print "CO.0003623 splice dict 2", ("CO.0003623" in splice_dict_2)	

	splice_lib.add_transcripts_to_event_dict(
		args.event_ioe_2, 
		splice_dict_2)



	transcript_dict_1 = splice_lib.generate_standard_transcript_dict(
											transcript_gtf_1)
	transcript_dict_2 = splice_lib.generate_standard_transcript_dict(
											transcript_gtf_2)



	splice_lib.sort_transcript_dict_exons(
						transcript_dict_1)
	splice_lib.add_junctions_to_transcript_dict(
						transcript_dict_1)
	splice_lib.sort_transcript_dict_exons(
						transcript_dict_2)
	splice_lib.add_junctions_to_transcript_dict(
						transcript_dict_2)




	transcript_match_dict = import_transcript_matching_dict(
								transcript_table,
								transcript_dict_1)



	species_1_transcript_centric_included_forms = make_transcript_event_dict(
														splice_dict_1, 
														"included")
	species_1_transcript_centric_excluded_forms = make_transcript_event_dict(
														splice_dict_1, 
														"excluded")
	species_2_transcript_centric_included_forms = make_transcript_event_dict(
														splice_dict_2, 
														"included")
	species_2_transcript_centric_excluded_forms = make_transcript_event_dict(
														splice_dict_2, 
														"excluded")

	matches = find_matches(
		splice_dict_1, 
		splice_dict_2, 
		transcript_match_dict, 
		transcript_dict_1, 
		transcript_dict_2, 
		species_1_transcript_centric_included_forms, 
		species_1_transcript_centric_excluded_forms, 
		species_2_transcript_centric_included_forms, 
		species_2_transcript_centric_excluded_forms)


	with open(
			outdir + 
			"/annotated_event_matches.tsv", 
			'w') as file:

		for event, event_match in matches.iteritems():

			file.write(event + 
					   "\t" + 
					   event_match["match"] + 
					   "\t" + 
					   event_match["support"] +
					   "\n")
			

if __name__ == '__main__':

	main()