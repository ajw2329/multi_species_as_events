'''
Takes as input find_switch_events.py output
from multiple runs and outputs a new file
with the consensus result.  For instance,
an event which is 'always' nmd in one input file
and 'sometimes' nmd in another becomes 'sometimes' 
overall.  Note that the input files must contain
all of the same event ids.
'''

import argparse
import sys
import subprocess



def process_status_form_info(
		status,
		form):

	'''
	Jointly process status and form to ensure that
	the resulting consensus makes sense for both.

	For instance, an nmd event could be "always"
	in all inputs but have different forms 
	in different inputs.  If this is the case, the
	"always" status no longer makes sense and must
	be changed to "sometimes"

	Note that currently 'coding_status' actually
	violates the consensus that 'always', 'sometimes',
	'never' refer to switching events and thus should not
	be checked using this function.
	'''

	consensus_status = get_consensus_category(status)

	consensus_form = get_consensus_category(form)


	if consensus_status == "always" and consensus_form == "both":

		consensus_status = "sometimes"





def get_consensus_category(category_set):

	if len(category_set) == 1:

		consensus_category = category_set.pop()


	elif not {"always","sometimes","never"}.isdisjoint(category_set):

		if {"never","NA"} == category_set:

			consensus_category = "never"

		else:

			consensus_category = "sometimes"


	elif not {"included","excluded","both"}.isdisjoint(category_set):

		if {"neither","NA"} == category_set:

			consensus_category = "neither"

		else:

			consensus_category = "both"

	else:

		sys.exit("Some category combination not properly accounted for by get_consensus_category()")

	return consensus_category




def parse_arguments(input_args):

	pass


def main(args):

	pass


if __name__ == "__main__":

	output_args = parse_arguments(sys.argv[1:])
	main(output_args)