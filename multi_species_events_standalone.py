import copy
import splice_lib
import argparse
import glob
import subprocess
import sys
import re
import cPickle as pickle
import gen_methods as gm
import infer_pairwise_events
import assign_events_to_genes


##TODO: Haircuts.  Lack of them is probably causing events to get purged due to slightly misplaced outer nodes
##TODO: Add some functionality to check whether exons from same event share strand/chrom when importing GTF 
##TODO: liftOver/minimap2 comparison strategy.  e.g. #events paired up that were identified in both species (particularly when this is semi-verifiable e.g. by using CAT gtfs), 
##TODO: figure out better way of dealing with non-reciprocal (but still multi-way) mapping issue (currently a warning is output)



def generate_minimap2_pre_combination_event_dict(event_gtf_filename):
    '''
        Generates a dictionary of events indexed by event ID
        {"SE.0000056":
            "event_type": "SE",
            "included_forms": {
                "A": {
                    "chrom": "chr1",
                    "strand": "+",
                    "exons": [[100, 200], [300, 400],[500,600]]
                },
            "excluded_forms": {
                "A": {
                    "chrom": "chr1",
                    "strand": "+",
                    "exons": [[100, 200],[500, 600]]
                },
                "B": {
                    "chrom": "chr5",
                    "strand": "-",
                    "exons": [[300, 400]]
                }
            }
        }

    '''

    event_types = ["A3", "A5", "AF", "AL", "RI", "MX", "SE", "MS", "MF", "ML", "CF", "CL", "CO", "AT", "AP", "MR"]

    standard_event_dict = {}

    pre_combination_event_dict = {}

    if event_gtf_filename.endswith(".gz"):

        gtf_file = gzip.open(event_gtf_filename, 'rb')

    else:

        gtf_file = open(event_gtf_filename, 'r')

    for line in gtf_file:

        entry = line.split()

        if entry[2] == "exon":

            start = int(entry[3])
            end = int(entry[4])

            chrom = re.sub("_","&",entry[0].strip())
            strand = entry[6].strip()

            transcript_id_entry = re.findall('transcript_id\s\"[^;\"]+\";', line)

            if len(transcript_id_entry) == 0:

                sys.exit("Bad gtf file - 'transcript_id' not found in exon entry")

            elif len(transcript_id_entry) > 1:

                print "Bad transcript_id - setting transcript_id to 'bad_transcript_id'"

                transcript_id = "bad_transcript_id"

            else:

                    transcript_id = re.sub('[";]', '', transcript_id_entry[0].strip().split()[1])    

 

            event_id = "_".join(transcript_id.split("_")[0:-1])
            base_event_id = event_id.split("&")[0] ## extract event ID without alpha numeric suffix
            version = event_id.split("&")[1]

            form = transcript_id.split("_")[-1]

            for i in event_types:

                if i in event_id:

                    event_type = i
                    break


            if base_event_id not in pre_combination_event_dict:

                pre_combination_event_dict[base_event_id] = {
                    "event_type": event_type,
                    "included_forms": {},
                    "excluded_forms": {}
                }

            if form == "included":

                pre_combination_event_dict[base_event_id]["included_forms"].setdefault(version, {"chrom": chrom, "strand": strand}).setdefault("exons", []).append([start, end])

            elif form == "excluded":

                pre_combination_event_dict[base_event_id]["excluded_forms"].setdefault(version, {"chrom": chrom, "strand": strand}).setdefault("exons", []).append([start, end])

            else:

                sys.exit("Non-standard event form found in generate_minimap2_pre_combination_event_dict of multi_species_events.py.  Exiting . . . ")


    gtf_file.close()

    return pre_combination_event_dict


def resolve_minimap_event_combinations(pre_combination_event_dict):


    ######Seems to be a major issue with this function - at the end of the day the choice of combinations is arbitrary.  No more than one (essentially random) combination is retained.
    ######Ideally any events with more than one possibility should be logged


    event_type_length_dict = {"MS": {"inc": 4, "exc": 2},
                              "SE": {"inc": 3, "exc": 2},
                              "A3": {"inc": 2, "exc": 2},
                              "A5": {"inc": 2, "exc": 2},
                              "AF": {"inc": 2, "exc": 2},
                              "AL": {"inc": 2, "exc": 2},
                              "RI": {"inc": 1, "exc": 2},
                              "MR": {"inc": 1, "exc": 3},
                              "CO": {"inc": 2, "exc": 2},
                              "MF": {"inc": 2, "exc": 2},
                              "ML": {"inc": 2, "exc": 2},
                              "CF": {"inc": 2, "exc": 2},
                              "CL": {"inc": 2, "exc": 2},
                              "MX": {"inc": 3, "exc": 3}}

    minimap2_standard_event_dict = {}

    for event in pre_combination_event_dict:

        event_type = pre_combination_event_dict[event]["event_type"]

        if len(pre_combination_event_dict[event]["included_forms"]) == 0:

            if len(pre_combination_event_dict[event]["excluded_forms"]) == 0: ######Is this a possible condition at all?  If not, remove it!!!

                continue

            else:

                if event_type is "RI":

                    for exc_version in pre_combination_event_dict[event]["excluded_forms"]:

                        exc_chrom = pre_combination_event_dict[event]["excluded_forms"][exc_version]["chrom"]
                        exc_strand = pre_combination_event_dict[event]["excluded_forms"][exc_version]["strand"]


                        exc_exons = pre_combination_event_dict[event]["excluded_forms"][exc_version]["exons"]

                        if len(exc_exons) == 2:

                            inc_exons = [[pre_combination_event_dict[event]["excluded_forms"][exc_version]["exons"][0][0], pre_combination_event_dict[event]["excluded_forms"][exc_version]["exons"][1][1]]]

                            ######This is actually eliminating possibilities.  By using "event" as the indexing variable, any subsequent excluded form versions will overwrite this one.

                            minimap2_standard_event_dict[event] = {
                                "event_type": event_type,
                                "included_exons": copy.deepcopy(inc_exons),
                                "excluded_exons": copy.deepcopy(exc_exons),
                                "chrom": exc_chrom,
                                "strand": exc_strand
                            }


        elif len(pre_combination_event_dict[event]["excluded_forms"]) == 0:

            if event_type in ["SE", "MS"]:


                for inc_version in pre_combination_event_dict[event]["included_forms"]:

                    inc_chrom = pre_combination_event_dict[event]["included_forms"][inc_version]["chrom"]
                    inc_strand = pre_combination_event_dict[event]["included_forms"][inc_version]["strand"]

                    inc_exons = pre_combination_event_dict[event]["included_forms"][inc_version]["exons"]

                    if (event_type is "SE" and len(inc_exons) == 3) or (event_type is "MS" and len(inc_exons) > 3):

                        exc_exons = [[pre_combination_event_dict[event]["included_forms"][inc_version]["exons"][0][0], pre_combination_event_dict[event]["included_forms"][inc_version]["exons"][0][1]], [pre_combination_event_dict[event]["included_forms"][inc_version]["exons"][-1][0], pre_combination_event_dict[event]["included_forms"][inc_version]["exons"][-1][1]]]

                        ######Same situation as noted on line 183 - using event as index will cause subsequent included forms to overwrite this one.

                        minimap2_standard_event_dict[event] = {
                            "event_type": event_type,
                            "included_exons": copy.deepcopy(inc_exons),
                            "excluded_exons": copy.deepcopy(exc_exons),
                            "chrom": inc_chrom,
                            "strand": inc_strand
                        }


        else:

            for inc_version in pre_combination_event_dict[event]["included_forms"]:

                inc_chrom = pre_combination_event_dict[event]["included_forms"][inc_version]["chrom"]
                inc_strand = pre_combination_event_dict[event]["included_forms"][inc_version]["strand"]

                inc_exons = pre_combination_event_dict[event]["included_forms"][inc_version]["exons"]

                for exc_version in pre_combination_event_dict[event]["excluded_forms"]:

                    exc_chrom = pre_combination_event_dict[event]["excluded_forms"][exc_version]["chrom"]
                    exc_strand = pre_combination_event_dict[event]["excluded_forms"][exc_version]["strand"]

                    exc_exons = pre_combination_event_dict[event]["excluded_forms"][exc_version]["exons"]

                    def add_to_dict():

                        minimap2_standard_event_dict[event] = {
                                "event_type": event_type,
                                "included_exons": copy.deepcopy(inc_exons),
                                "excluded_exons": copy.deepcopy(exc_exons),
                                "chrom": inc_chrom,
                                "strand": inc_strand
                            }

                    if inc_chrom == exc_chrom and inc_strand == exc_strand:

                        if event_type == "MS":

                            if len(inc_exons) >= event_type_length_dict[event_type]["inc"] and len(exc_exons) == event_type_length_dict[event_type]["exc"]:

                                add_to_dict()

                        elif event_type == "MR":

                            if len(inc_exons) == event_type_length_dict[event_type]["inc"] and len(exc_exons) >= event_type_length_dict[event_type]["exc"]:

                                add_to_dict()

                        elif event_type in ["SE", "RI", "A3", "A5", "AF", "AL", "MX"]:

                            if len(inc_exons) == event_type_length_dict[event_type]["inc"] and len(exc_exons) == event_type_length_dict[event_type]["exc"]:

                                add_to_dict()

                        elif event_type in ["CO", "CF", "CL"]:

                            if len(inc_exons) >= 2 and len(exc_exons) >= 2:

                                add_to_dict()

                        elif event_type in ["MF", "ML"]:

                            if (len(inc_exons) >= 3 and len(exc_exons) >= 2) or (len(inc_exons) >= 2 and len(exc_exons) >= 3):

                                add_to_dict()


    return minimap2_standard_event_dict



def deduplicated_minimap2_sam(sam_outfile_path, newsam_path):

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    new_samfile = open(newsam_path, 'w')

    used_IDs = {}

    with open(sam_outfile_path, 'r') as file:

        for line in file:

            if not line.startswith("@"):

                entry = line.split()

                ID = entry[0].strip() 


                if ID in used_IDs:

                    used_IDs[ID] += 1

                    new_ID = ID.split("_")[0] + "&" + alphabet[used_IDs[ID]] + "_" + ID.split("_")[1]
                    entry[0] = new_ID

                else:

                    used_IDs[ID] = 0

                    new_ID = ID.split("_")[0] + "&" + alphabet[used_IDs[ID]] + "_" + ID.split("_")[1]
                    entry[0] = new_ID

                new_samfile.write("\t".join(entry) + "\n")

            else:

                new_samfile.write(line.strip() + "\n")

    new_samfile.close()



def get_cross_species_events(from_species, to_species, species_dict, num_threads, minimap2_path = "minimap2", samtools_path = "samtools", bedtools_path = "bedtools", bedToGenePred_path = "bedToGenePred", genePredToGtf_path = "genePredToGtf"):

    #/hive/users/anjowall/minimap2/minimap2 -a -G 500000 -I 100G -t 10 -u f -x splice /hive/users/anjowall/genomes/rheMac8/STAR_genome_rheMac8/rheMac8.fa splice_lib_events.fa > minimap2_test/hu_to_rh_test.sam

    sam_outfile_path = species_dict[to_species]["minimap_outdir"] + "/" + to_species + "_" + "minimapped_from_" + from_species + "_splice_events.sam"

    newsam_path = sam_outfile_path.split(".sam")[0] + ".dedup_id.sam"


    try: ##from Paul Finn http://www.pfinn.net/python-check-if-file-exists.html
        with open(newsam_path[:-3] + "gtf") as file:

            print newsam_path[:-3] + "gtf exists - skipping minimap run."

    except IOError:

        try:

            subprocess.check_output(minimap2_path + " -a -G 300000 --secondary=no -t " + num_threads + " -u f -x splice " + species_dict[to_species]["genome_fasta"] + " " + species_dict[from_species]["event_fasta_path"] + " > " + sam_outfile_path, shell = True, stderr = subprocess.STDOUT)

            deduplicated_minimap2_sam(sam_outfile_path, newsam_path)

            subprocess.check_output(samtools_path + " view -hb " + newsam_path + " > " + newsam_path[:-3] + "bam", shell = True, stderr = subprocess.STDOUT)
            subprocess.check_output(samtools_path + " sort " + newsam_path[:-3] + "bam" + " > " + newsam_path[:-3] + "sorted.bam", shell = True, stderr = subprocess.STDOUT)
            subprocess.check_output(samtools_path + " index " + newsam_path[:-3] + "sorted.bam", shell = True, stderr = subprocess.STDOUT)
            subprocess.check_output(bedtools_path + " bamtobed -bed12 -i " + newsam_path[:-3] + "sorted.bam" + " > " + newsam_path[:-3] + "bed12", shell = True, stderr = subprocess.STDOUT)
            subprocess.check_output(bedToGenePred_path + " " + newsam_path[:-3] + "bed12 stdout | " + genePredToGtf_path + " file stdin stdout | awk '$3==\"exon\"' > " + newsam_path[:-3] + "gtf", shell = True, stderr = subprocess.STDOUT)

        ##import as splice_lib event dict --- this needs something special to parse modified event id


        ###Need some code to deal with orphaned forms - e.g. orphaned SE included form actually serves as entire form etc etc (evaluate_orphans function)

        ###Also run some complete event dict code

        except subprocess.CalledProcessError as e:
            print e.output
            sys.exit("minimap2 failed in 'get_cross_species_de_novo_events' of multi_species_events_standalone.py. Exiting . . . ")

    precombination_event_dict = generate_minimap2_pre_combination_event_dict(newsam_path[:-3] + "gtf")
    minimap2_standard_event_dict = resolve_minimap_event_combinations(precombination_event_dict)

    splice_lib.complete_event_dict(minimap2_standard_event_dict)

    return minimap2_standard_event_dict





def check_event_type(event_dict_entry, event):

    '''
    Reconstructs input for and calls 'classify_event' from infer_pairwise_events.py.  Return value is a tuple consisting of event_type
    
    '''

    ######Try to change this in order to account for the possibility that events just need a haircut
    ######Write haircut method to perform all the required checks

    event_misclassified = False

    pre_common = []
    post_common = []
    tx1_unique = []
    tx2_unique = []
    form_1_exons = copy.deepcopy(event_dict_entry["included_exons"])
    form_2_exons = copy.deepcopy(event_dict_entry["excluded_exons"])
    strand = event_dict_entry["strand"]
    tx1 = "tx1"
    tx2 = "tx2"

    form_1_exons = [map(int, i) for i in form_1_exons]
    form_2_exons = [map(int, i) for i in form_2_exons]

    def convert_to_node_list(exons):

        nodes = []

        for exon in exons:

            nodes.append(str(exon[0]) + "_left")
            nodes.append(str(exon[1]) + "_right")

        return nodes

    flattened_form_1 = convert_to_node_list(form_1_exons)
    flattened_form_2 = convert_to_node_list(form_2_exons)

    tx1_unique = [i for i in flattened_form_1 if i not in flattened_form_2]
    tx2_unique = [i for i in flattened_form_2 if i not in flattened_form_1]

    if len(tx1_unique) == 0 and len(tx2_unique) == 0:

        return False

    while flattened_form_1[0] == flattened_form_2[0]:

        pre_common.append(flattened_form_1.pop(0))
        flattened_form_2.pop(0)

        if len(flattened_form_1) == 0 or len(flattened_form_2) == 0:

            break

    if len(flattened_form_1) > 1:

        while flattened_form_1[0] in tx1_unique:

            flattened_form_1.pop(0)
            if len(flattened_form_1) == 0:
                break

    if len(flattened_form_2) > 0:
    
        while flattened_form_2[0] in tx2_unique:

            flattened_form_2.pop(0)
            if len(flattened_form_2) == 0:
                break

    if len(flattened_form_1) != 0 and len(flattened_form_2) != 0:

        while flattened_form_1[0] == flattened_form_2[0]:

            post_common.append(flattened_form_1.pop(0))
            flattened_form_2.pop(0)

            if len(flattened_form_1) == 0 or len(flattened_form_2) == 0:

                break

    initial_event_type = event_dict_entry["event_type"]


    if (initial_event_type in ["AF","MF","CF","UF"] and strand == "+") or initial_event_type in ["AL","ML","CL","UL"] and strand == "-":

        if event_dict_entry["included_exons"][0][0] != event_dict_entry["excluded_exons"][0][0] and event_dict_entry["included_exons"][-1][-1] == event_dict_entry["excluded_exons"][-1][-1]:

            form_1_exons.insert(0, "universal_left")
            form_2_exons.insert(0, "universal_left")
            pre_common.append("universal_left")

        else:

            print event, "Pre-comprehensive misclassification tag"
            event_misclassified = True ### replace this with something more helpful


    elif (initial_event_type in ["AF","MF","CF","UF"] and strand == "-") or initial_event_type in ["AL","ML","CL","UL"] and strand == "+":

        if event_dict_entry["included_exons"][-1][-1] != event_dict_entry["excluded_exons"][-1][-1] and event_dict_entry["included_exons"][0][0] == event_dict_entry["excluded_exons"][0][0]:

            form_1_exons.append("universal_right")
            form_2_exons.append("universal_right")
            post_common.append("universal_right")


        else:

            event_misclassified = True ### replace this with something more helpful
            print event, "Pre-comprehensive misclassification tag"

    else:

        if not (event_dict_entry["included_exons"][0][0] == event_dict_entry["excluded_exons"][0][0] and event_dict_entry["included_exons"][-1][-1] == event_dict_entry["excluded_exons"][-1][-1]):

            event_misclassified = True ### replace this with something more helpful
            print event, "Pre-comprehensive misclassification tag"

    if not event_misclassified:

        class_result = infer_pairwise_events.classify_event(pre_common, post_common, tx1_unique, tx2_unique, form_1_exons, form_2_exons, strand, tx1, tx2)

        print event, class_result

        same_classification = initial_event_type == class_result["event_type"]

        return same_classification

    else:

        return False


def get_junction_key(event_dict_entry):

    flat_joined_included_jl = "_".join([str(i) for j in event_dict_entry["included_junctions"] for i in j])
    flat_joined_excluded_jl = "_".join([str(i) for j in event_dict_entry["excluded_junctions"] for i in j])

    flat_joined_included_eij = "_".join(map(str, event_dict_entry["included_ei_junctions"]))
    flat_joined_excluded_eij = "_".join(map(str, event_dict_entry["excluded_ei_junctions"]))

    chrom = event_dict_entry["chrom"]
    strand = event_dict_entry["strand"]

    key = chrom + "_" + flat_joined_included_jl + "|" + flat_joined_included_eij + ":" +  flat_joined_excluded_jl + "|" + flat_joined_excluded_eij + "_" + strand



def merged_event_dicts_by_species(event_dicts, cross_species_event_dicts, species_in_order, outdir):

    possible_duplicate_log = open(outdir + "/possible_duplicate_log.txt", 'w')

    merged_events = {}
    renamed_merged_events = {}

    for species in event_dicts:

        merged_events[species] = {}

        for event in event_dicts[species]:

            junction_key = get_junction_key(event_dicts[species][event])

            merged_events[species][junction_key] = copy.deepcopy(event_dicts[species][event])

            merged_events[species][junction_key]["native"] = True
            merged_events[species][junction_key]["native_id"] = event
            merged_events[species][junction_key]["cross_species_keys"] = {}
            merged_events[species][junction_key]["cross_species_keys"][species] = {event: junction_key}


        for other_species in cross_species_event_dicts[species]:

            for event in cross_species_event_dicts[species][other_species]:

                junction_key = get_junction_key(cross_species_event_dicts[species][other_species][event])


                if event in event_dicts[other_species]:

                    other_junction_key = get_junction_key(event_dicts[other_species][event])

                else:

                    print >> sys.stderr, event, "Not found in event dict of", other_species
                    other_junction_key = None

                if junction_key in merged_events[species]:

                    ###Must reference the junction key from the source species !!!!!!

                    if event in event_dicts[other_species]:

                        if other_junction_key:

                            merged_events[species][junction_key]["cross_species_keys"][other_species] = {event: other_junction_key}

                else:

                    merged_events[species][junction_key] = copy.deepcopy(cross_species_event_dicts[species][other_species][event])
                    merged_events[species][junction_key]["native"] = False
                    merged_events[species][junction_key]["native_id"] = None

                    if other_junction_key:
                        merged_events[species][junction_key]["cross_species_keys"] = {other_species: {event: other_junction_key}, species: {"placeholder": junction_key}}
                    else:
                        merged_events[species][junction_key]["cross_species_keys"] = {species: {"placeholder": junction_key}}

    
    i = 2
    while i > 0:

        for species in merged_events:

            #print species, "1"

            for junction_key in merged_events[species]:

                #print junction_key, "2"
                master_cross_species_keys = {}

                if species in merged_events[species][junction_key]["cross_species_keys"]:

                    event = merged_events[species][junction_key]["cross_species_keys"][species].keys()[0]

                else:

                    event = "placeholder"
                    merged_events[species][junction_key]["cross_species_keys"][species] = {event: junction_key}

                for other_species in merged_events[species][junction_key]["cross_species_keys"]:

                    #print other_species, "3"

                    other_id = merged_events[species][junction_key]["cross_species_keys"][other_species].keys()[0]
                    other_junction_key =  merged_events[species][junction_key]["cross_species_keys"][other_species][other_id]

                    #print other_id, "4"
                    #print other_junction_key, "5"

                    if other_junction_key in merged_events[other_species]:

                        #print "yes", "6"

                        master_cross_species_keys.update(merged_events[other_species][other_junction_key]["cross_species_keys"])
                        merged_events[other_species][other_junction_key]["cross_species_keys"].update(merged_events[species][junction_key]["cross_species_keys"])

                    else:
                        possible_duplicate_log.write("Please investigate event with putative junction key " + other_junction_key + " in species " + other_species + " which putatively connects to event with junction key " + junction_key + " in species " + species + "\n")

                else:
                    merged_events[species][junction_key]["cross_species_keys"].update(master_cross_species_keys)

        i -= 1

    ####Name event by human ID whenever possible.  Species-specific events should be denoted as such in ID.  e.g. cSE.0000457 for chimpanzee-specific SE, or cmoSE.0000301 for SE event common to chimpanzee, macaque, orangutan

    #pickle.dump(merged_events, open(outdir + "/merged_events_pre_rename_standard_event_dict.pkl", 'wb'), pickle.HIGHEST_PROTOCOL)

    event_counter = 1

    all_break = False

    for species in species_in_order:

        if species in merged_events:

            while len(merged_events[species]) > 0:
                
                if all_break:
                    
                    break

                working_event_key = merged_events[species].keys()[0]

                all_species_keys = {}

                species_in_id = []

                for ref_species in species_in_order:

                    if ref_species in merged_events[species][working_event_key]["cross_species_keys"]:

                        #ref_species_id = merged_events[species][junction_key]["cross_species_keys"][ref_species].keys()[0]
                        ref_species_junction_key = merged_events[species][working_event_key]["cross_species_keys"][ref_species][merged_events[species][working_event_key]["cross_species_keys"][ref_species].keys()[0]]

                        all_species_keys[ref_species] = ref_species_junction_key

                        species_in_id.append(ref_species)

                new_id_primary = merged_events[species][working_event_key]["event_type"] + "." + str(event_counter).zfill(7)
                full_new_id = "_".join(species_in_id) + "|" + new_id_primary

                for ref_species in all_species_keys:

                    if all_species_keys[ref_species] in merged_events[ref_species]:

                        #print full_new_id
                        #print ref_species, all_species_keys[ref_species]
                        renamed_merged_events.setdefault(ref_species, {}).setdefault(full_new_id, copy.deepcopy(merged_events[ref_species][all_species_keys[ref_species]]))

                        del merged_events[ref_species][all_species_keys[ref_species]]

                    else:

                        possible_duplicate_log.write("Please investigate event with putative junction key " + all_species_keys[ref_species] + " in species " + ref_species + " which putatively connects to event with junction key " + "\n")
                        possible_duplicate_log.write(full_new_id + " has possible duplicate in " + species + " with key " + working_event_key + "\n")
                        
                else:
                    
                    if working_event_key in merged_events[species]:
                        
                        del merged_events[species][working_event_key]
     

                event_counter += 1

    possible_duplicate_log.close()

    return renamed_merged_events


def complete_merged_dict(merged_event_dict, species_dict, add_gene_symbols, transcript_gtfs):

    ##add junctions, junction counts slots, remove duplicates or otherwise problematic events, clean up outermost exon edges


    for species in merged_event_dict:

        splice_lib.complete_event_dict(merged_event_dict[species])

        if add_gene_symbols:

            if transcript_gtfs:

                matches_dict, genes_dict, transcript_dict = assign_events_to_genes.main(["--transcript_gtf", species_dict[species]["merged_stringtie_gtf_path"], "--gene_gtf", species_dict[species]["annotated_transcript_gtf"], "--prefix", "placeholder", "--suppress_output"], merged_event_dict[species])

            else:

                matches_dict, genes_dict, transcript_dict = assign_events_to_genes.main(["--gene_gtf", species_dict[species]["annotated_transcript_gtf"], "--prefix", "placeholder", "--suppress_output"], merged_event_dict[species])

            for gene in matches_dict:

                for event in matches_dict[gene]:

                    merged_event_dict[species][event].setdefault("gene", []).append(gene)




def check_gene_equivalence(merged_event_dict, outdir, species_in_order):

    '''
        Check that an event maps to the equivalent gene (using CAT annotations or similar) in all species.  If it doesn't, write the information to an output file so it can be investigated and/or discarded in downstream analysis
    '''

    events_genes = {}

    multiple_genes = open(outdir + "/events_with_contrary_genes.tsv", 'w')

    multiple_genes.write("\t".join(["event_id"] + species_in_order + ["consensus_overlap"]) + "\n")


    for species in merged_event_dict:

        for event in merged_event_dict[species]:

            if "gene" in merged_event_dict[species][event]:

                events_genes.setdefault(event, {}).update({species: ",".join(merged_event_dict[species][event]["gene"])})

    for event in events_genes:

        all_genes = []

        for species in species_in_order:

            if species in events_genes[event]:

                all_genes.append(events_genes[event][species])


        if len(set(all_genes)) > 1:

            consensus_overlap = False

            if len(set.intersection(*[set(i.split(",")) for i in all_genes])) > 0:

                consensus_overlap = True

            genes_to_write = []

            for species in species_in_order:

                if species in events_genes[event]:

                    genes_to_write.append(events_genes[event][species])

                else:

                    genes_to_write.append("NA")

            multiple_genes.write("\t".join([event] + genes_to_write + [str(consensus_overlap)]) + "\n")

    multiple_genes.close()






def run_liftOver(path_to_chain_file, input_gtf_path, outdir, from_species, to_species):

    path_to_successful_lifts = outdir + "/" + from_species + "_to_" + to_species + "_lift.gtf"

    try: ##from Paul Finn http://www.pfinn.net/python-check-if-file-exists.html
        with open(path_to_successful_lifts) as file:

            print path_to_successful_lifts + " gtf exists - skipping liftOver run."

    except IOError:

        try:

            cmd_str = "liftOver -gff " + input_gtf_path + " " + path_to_chain_file + " " + path_to_successful_lifts + " " + outdir + "/" + from_species + "_to_" + to_species + "_lift_fail.gtf"

            gm.run_bash_cmd(cmd_str.split())

        except subprocess.CalledProcessError as e:
            print e.output
            sys.exit("liftOver failed in 'run_liftOver' of multi_species_events_standalone.py. Exiting . . . ")

    return path_to_successful_lifts


def same_exon_counts(event_dict_entry_list):
    '''
        takes as input a list of event dict entries and checks whether they all have the same number of included and excluded form exons.  Returns boolean.
    '''

    included_exon_counts = [len(i["included_exons"]) for i in event_dict_entry_list]
    excluded_exon_counts = [len(i["excluded_exons"]) for i in event_dict_entry_list]

    if len(set(included_exon_counts)) > 1 or len(set(excluded_exon_counts)) > 1:

        return False

    return True




def check_mapped_event_integrity(mapped_species, source_species, mapped_events_dict, source_events_dict, mapping_failure_log_file):
    '''
        Ensures that mapped events have same exon counts for both source and mapped features, as well as the same event type

        Assumes that mapped_events_dict and source_events_dict have same keys
    '''

    event_issues = {}

    for event in mapped_events_dict:

        write_to_log = False

        if not same_exon_counts([mapped_events_dict[event], source_events_dict[event]]):

            event_issues.setdefault(event, set()).add("exon_count")
            write_to_log = True


        if not (len(mapped_events_dict[event]["included_exons"]) == 0 or len(mapped_events_dict[event]["excluded_exons"]) == 0):

            print mapped_events_dict[event]["included_exons"]
            print mapped_events_dict[event]["excluded_exons"]

            same_classification = check_event_type(mapped_events_dict[event], event)

            if not same_classification:

                event_issues.setdefault(event, set()).add("post_map_class_error")
                write_to_log = True

        else:

            same_classification = False
            event_issues.setdefault(event, set()).add("post_map_class_error")
            write_to_log = True

        if write_to_log:

            mapping_failure_log_file.write("\t".join([mapped_species + "," + source_species, event, ",".join(event_issues[event])]) + "\n")


    for event in event_issues:

        del mapped_events_dict[event]


def final_original_name(merged_event_dict, species_list, outdir):

    #merged_events[species][junction_key]["cross_species_keys"][species]


    species_comparisons = []

    for i, species_1 in enumerate(species_list[:-1]):

        for j, species_2 in enumerate(species_list[i+1:]):

            species_comparisons.append([species_1,species_2])

    entry_set = set()

    for species in merged_event_dict:

        for event in merged_event_dict[species]:

            for comparison in species_comparisons:

                if comparison[0] in merged_event_dict[species][event]["cross_species_keys"] and comparison[1] in merged_event_dict[species][event]["cross_species_keys"]:


                    if len(merged_event_dict[species][event]["cross_species_keys"][comparison[0]].keys()) > 1 or len(merged_event_dict[species][event]["cross_species_keys"][comparison[1]].keys()) > 1:

                        print "WARNING: event", event, "has more than one cross-species key for species", comparison[0], "or", comparison[1]

                    entry_string = "\t".join([event, merged_event_dict[species][event]["cross_species_keys"][comparison[0]].keys()[0], merged_event_dict[species][event]["cross_species_keys"][comparison[1]].keys()[0], comparison[0], comparison[1]])

                    entry_set.add(entry_string)



    with open(outdir + "/final_original_name.tsv", 'w') as file:

        file.write("\t".join(["final_event_id", "original_event_id_species_1", "original_event_id_species_2", "species_1", "species_2"]) + "\n")

        for entry in entry_set:

            file.write(entry + "\n")


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--event_file_prefixes", type = str, help = "Comma-separated list of prefixes to event gtfs, ioes, and fasta files.  e.g. human_splice_lib_events,chimpanzee_splice_lib_events which will allow the program to find human_splice_lib_events.gtf, human_splice_lib_events.fa, chimpanzee_splice_lib_events.gtf and chimpanzee_splice_lib_events.fa from input directory", required = True)
    parser.add_argument("--species_list", type = str, help = "Comma-separated list of species names.  Expected to be in same order as corresponding prefixes supplied in --event_file_prefixes argument", required = True)
    parser.add_argument("--indir", type = str, help = "Path to input directory in which event gtf and fasta files can be found (.fa extension expected for fasta files) with prefixes matching those supplied to the --event_file_prefixes argument", required = True)
    parser.add_argument("--genome_fasta_paths", type = str, help = "Comma-separated paths to genome fasta files.  In same order as '--species_list' argument")
    parser.add_argument("--chain_paths", type = str, help = "Comma-separated paths to chain files.  For each pair of species, a chain file should be supplied for both directions.  So, for n-species, n(n-1) chain files should be provided.  These should be sorted first by the 'from' species, and next by the 'to' species, with each sorting in same order as '--species_list' argument.  E.g. for three species such as human,chimp,orang, the chain paths should be ordered as follows: /path/to/human_to_chimp,/path/to/human_to_orang,/path/to/chimp_to_human,/path/to/chimp_to_orang,/path/to/orang_to_human,/path/to/orang_to_chimp")
    parser.add_argument("--minimap_outdir", type = str, help = "Path to directory for minimap output.  Providing this argument also specifies that you want to run minimap2.")
    parser.add_argument("--liftover_outdir", type = str, help = "Path to directory for liftover output.  Providing this argument also specifies that you want to run liftOver")
    parser.add_argument("--samtools_path", type = str, help = "Path to samtools executable - default = 'samtools'", default = "samtools")
    parser.add_argument("--bedtools_path", type = str, help = "Path to bedtools executable - default = 'bedtools'", default = "bedtools")
    parser.add_argument("--bedToGenePred_path", type = str, help = "Path to bedToGenePred executable - default = 'bedToGenePred'", default = "bedToGenePred")
    parser.add_argument("--minimap2_path", type = str, help = "Path to minimap2 executable - default = 'minimap2'", default = "minimap2")
    parser.add_argument("--num_threads", type = str, help = "Number of threads for minimap2 to use - default = 1", default = "1")
    parser.add_argument("--genePredToGtf_path", type = str, help = "Path to genePredToGtf executable - default = 'genePredToGtf'", default = "genePredToGtf")
    parser.add_argument("--liftOver_path", type = str, help = "Path to liftOver executable - default = 'liftOver'", default = "liftOver")
    parser.add_argument("--add_gene_symbols", action = "store_true", help = "If set, the genes to which the events belong will be included.  Requires --gene_gtfs and can make use of --transcript_gtfs")
    parser.add_argument("--gene_gtfs", type = str, help = "Comma-separated list of paths to annotation gtfs for each species, in the same order as --species_list.  Required (but only used for) --add_gene_symbols")
    parser.add_argument("--transcript_gtfs", type = str, help = "Comma-separated list of paths to transcript gtfs for each species in the same order as --species_list.  Can be made use of by --add_gene_symbols to increase the number of events for which symbols can be identified")

    args = parser.parse_args()

    if args.add_gene_symbols and not args.gene_gtfs:

        sys.exit("Passage of --add_gene_symbols requires --gene_gtfs")

    elif args.add_gene_symbols:

        annotated_transcript_gtf_list = args.gene_gtfs.split(",")

        if args.transcript_gtfs:

            merged_stringtie_gtf_path_list = args.transcript_gtfs.split(",")

    ##require that either minimap_outdir OR liftover_outdir be set

    if args.minimap_outdir is None and args.liftover_outdir is None:

        sys.exit("Neither --minimap_outdir nor --liftover_outdir are set.  At least one must be set")


    event_dicts = {}
    cross_species_event_dicts = {}
    species_dict = {}

    species_list = args.species_list.split(",")
    prefixes = args.event_file_prefixes.split(",")


    ## Make output directory

    if args.minimap_outdir:

        if args.genome_fasta_paths is None:
            sys.exit("Use of minimap requires genome fasta files.  Pass comma-separated paths to --genome_fasta_paths")

        gm.run_bash_cmd(("mkdir -p " + args.minimap_outdir).split())
        genome_fasta_paths = args.genome_fasta_paths.split(",")


    if args.liftover_outdir:

        if args.chain_paths is None:
            sys.exit("Use of liftover requires chain files.  Pass comma-separated paths to --chain_paths")

        gm.run_bash_cmd(("mkdir -p " + args.liftover_outdir).split())
        chain_paths = args.chain_paths.split(",")
        ###place chain file paths in dict indexed first by from species, then by to species

        chain_paths_dict = {}


        for species in species_list:

            species_list_minus = [i for i in species_list if i != species]

            for other_species in species_list_minus:

                chain_paths_dict.setdefault(species, {}).update({other_species: chain_paths.pop(0)})



    for i, species in enumerate(species_list):

        species_dict[species] = {}

        if args.minimap_outdir:
            species_outdir = args.minimap_outdir + "/" + species + "_out/"
            species_dict[species]["minimap_outdir"] = species_outdir
            gm.run_bash_cmd(("mkdir -p " + species_outdir).split())
            species_dict[species]["genome_fasta"] = genome_fasta_paths[i]

        if args.liftover_outdir:
            species_outdir = args.liftover_outdir + "/" + species + "_out/"
            species_dict[species]["liftover_outdir"] = species_outdir
            gm.run_bash_cmd(("mkdir -p " + species_outdir).split())

        if args.add_gene_symbols:
            species_dict[species]["annotated_transcript_gtf"] = annotated_transcript_gtf_list[i]
            if args.transcript_gtfs:
                species_dict[species]["merged_stringtie_gtf_path"] = merged_stringtie_gtf_path_list[i]


        species_dict[species]["prefix"] = prefixes[i]

        species_dict[species]["native_event_gtf_path"] = args.indir + "/" + species_dict[species]["prefix"] + ".gtf"

        try: 
            with open(species_dict[species]["native_event_gtf_path"]) as file:
                pass
        except IOError:
            sys.exit(species_dict[species]["native_event_gtf_path"] + " not found. Event GTF required for all species.  Exiting . . . ")


        if args.minimap_outdir:
            species_dict[species]["event_fasta_path"] = args.indir + "/" + species_dict[species]["prefix"] + ".fa"
            try:
                with open(species_dict[species]["event_fasta_path"]) as file:
                    pass
            except IOError:
                sys.exit(species_dict[species]["native_event_gtf_path"] + " not found. Event fasta required for all species for minimap2-based approach.  Exiting . . . ")


        event_dicts[species] = splice_lib.generate_standard_event_dict(species_dict[species]["native_event_gtf_path"])
        splice_lib.complete_event_dict(event_dicts[species])


    if args.minimap_outdir:

        mapping_failure_log = open(args.minimap_outdir + "/mapping_failure_log.tsv", 'w')
        mapping_failure_log.write("\t".join(["mapped_species,source_species", "event_id", "event_issues"]) + "\n")


        for i,species in enumerate(species_list):

            for other_species in species_list:

                if other_species != species:

                    cross_species_event_dicts.setdefault(other_species, {species: get_cross_species_events(species, other_species, species_dict, args.num_threads, minimap2_path = args.minimap2_path, samtools_path = args.samtools_path, bedtools_path = args.bedtools_path, bedToGenePred_path = args.bedToGenePred_path, genePredToGtf_path = args.genePredToGtf_path)}).update({species: get_cross_species_events(species, other_species, species_dict, args.num_threads, minimap2_path = args.minimap2_path, samtools_path = args.samtools_path, bedtools_path = args.bedtools_path, bedToGenePred_path = args.bedToGenePred_path, genePredToGtf_path = args.genePredToGtf_path)})

                    ##check for non-zero exon lists
                    ##check for preservation of exon counts
                    ##check for preservation of event type

                    check_mapped_event_integrity(other_species, species, cross_species_event_dicts[other_species][species], event_dicts[species], mapping_failure_log)


                    #splice_lib.check_integrity(cross_species_event_dicts[other_species][species])
                    splice_lib.complete_event_dict(cross_species_event_dicts[other_species][species])

        mapping_failure_log.close()
                

        merged_event_dict = merged_event_dicts_by_species(event_dicts, cross_species_event_dicts, species_list, args.minimap_outdir)

        complete_merged_dict(merged_event_dict, species_dict, args.add_gene_symbols, args.transcript_gtfs)

        native_events = {}


        for species in merged_event_dict:

            for event in merged_event_dict[species]:

                if merged_event_dict[species][event]["native"]:

                    native_events.setdefault(event, []).append(species)




        with open(args.minimap_outdir + "/native_event_species.tsv", "w") as file:

            file.write("\t".join(["event_id", "species_native"]) + "\n")

            for event in native_events:

                file.write("\t".join([event, ",".join(sorted(native_events[event]))]) + "\n")

        check_gene_equivalence(merged_event_dict, args.minimap_outdir, species_list)

        for species in merged_event_dict:

            splice_lib.output_event_gtf(merged_event_dict[species], args.minimap_outdir, name= species + "_splice_lib_events")
            splice_lib.output_event_bedfile(merged_event_dict[species], args.minimap_outdir, name = species + "_splice_lib_events")
            splice_lib.output_event_gtf_with_transcripts(merged_event_dict[species], args.minimap_outdir, name = species + "_splice_lib_events_with_transcripts")
            splice_lib.output_miso_event_gff3(merged_event_dict[species], args.minimap_outdir, name = species + "_splice_lib_events")


        final_original_name(merged_event_dict, species_list, args.minimap_outdir)



    if args.liftover_outdir:

        mapping_failure_log = open(args.liftover_outdir + "/mapping_failure_log.tsv", 'w')
        mapping_failure_log.write("\t".join(["mapped_species,source_species", "event_id", "event_issues"]) + "\n")

        for i,species in enumerate(species_list):

            for other_species in species_list:

                if other_species != species:

                    cross_species_event_dicts.setdefault(other_species, {species: splice_lib.generate_standard_event_dict(run_liftOver(chain_paths_dict[species][other_species], species_dict[species]["native_event_gtf_path"], args.liftover_outdir, species, other_species))}).update({species: splice_lib.generate_standard_event_dict(run_liftOver(chain_paths_dict[species][other_species], species_dict[species]["native_event_gtf_path"], args.liftover_outdir, species, other_species))})

                    ##check for non-zero exon lists
                    ##check for preservation of exon counts
                    ##check for preservation of event type

                    check_mapped_event_integrity(other_species, species, cross_species_event_dicts[other_species][species], event_dicts[species], mapping_failure_log)

                    splice_lib.complete_event_dict(cross_species_event_dicts[other_species][species])

        mapping_failure_log.close()



        merged_event_dict = merged_event_dicts_by_species(event_dicts, cross_species_event_dicts, species_list, args.liftover_outdir)

        complete_merged_dict(merged_event_dict, species_dict, args.add_gene_symbols, args.transcript_gtfs)

        native_events = {}


        for species in merged_event_dict:

            for event in merged_event_dict[species]:

                if merged_event_dict[species][event]["native"]:

                    native_events.setdefault(event, []).append(species)



        with open(args.liftover_outdir + "/native_event_species.tsv", "w") as file:

            file.write("\t".join(["event_id", "species_native"]) + "\n")

            for event in native_events:

                file.write("\t".join([event, ",".join(sorted(native_events[event]))]) + "\n")

        check_gene_equivalence(merged_event_dict, args.liftover_outdir, species_list)

        for species in merged_event_dict:

            splice_lib.output_event_gtf(merged_event_dict[species], args.liftover_outdir, name= species + "_splice_lib_events")
            splice_lib.output_event_bedfile(merged_event_dict[species], args.liftover_outdir, name = species + "_splice_lib_events")
            splice_lib.output_event_gtf_with_transcripts(merged_event_dict[species], args.liftover_outdir, name = species + "_splice_lib_events_with_transcripts")
            splice_lib.output_miso_event_gff3(merged_event_dict[species], args.liftover_outdir, name = species + "_splice_lib_events")

        final_original_name(merged_event_dict, species_list, args.liftover_outdir)



if __name__ == '__main__':

    main()
