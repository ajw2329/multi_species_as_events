import copy
from splice_lib import splice_lib
import argparse
import glob
import subprocess
import sys
import re
import cPickle as pickle
import gen_methods as gm
from junctionCounts import infer_pairwise_events
import assign_events_to_genes


##TODO: Haircuts.  Lack of them is probably causing events to get purged due to slightly misplaced outer nodes
##TODO: Add some functionality to check whether exons from same event share strand/chrom when importing GTF 
##TODO: liftOver/minimap2 comparison strategy.  e.g. #events paired up that were identified in both species (particularly when this is semi-verifiable e.g. by using CAT gtfs), 
##TODO: figure out better way of dealing with non-reciprocal (but still multi-way) mapping issue (currently a warning is output)



#### look into human_macaque-MF.0106250 from latest minimap run



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

    event_types = ["A3", 
                   "A5", 
                   "AF", 
                   "AL", 
                   "RI", 
                   "MX", 
                   "SE",
                   "MS", 
                   "MF", 
                   "ML", 
                   "CF", 
                   "CL", 
                   "CO", 
                   "AT", 
                   "AP", 
                   "MR"]

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

            transcript_id_entry = re.findall(
                'transcript_id\s\"[^;\"]+\";', 
                line)

            if len(transcript_id_entry) == 0:

                sys.exit("Bad gtf file - 'transcript_id' " +
                         "not found in exon entry")

            elif len(transcript_id_entry) > 1:

                print ("Bad transcript_id - setting " + 
                       "transcript_id to 'bad_transcript_id'")

                transcript_id = "bad_transcript_id"

            else:

                transcript_id = re.sub(
                    '[";]', '', 
                    transcript_id_entry[0].strip().split()[1])    

 

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

                ((pre_combination_event_dict[base_event_id]["included_forms"]
                    ).setdefault(version, {"chrom": chrom, "strand": strand})
                        ).setdefault("exons", []).append([start, end])

            elif form == "excluded":

                ((pre_combination_event_dict[base_event_id]["excluded_forms"]
                    ).setdefault(version, {"chrom": chrom, 
                                           "strand": strand})
                        ).setdefault("exons", []).append([start, end])

            else:

                sys.exit("Non-standard event form found " + 
                         "in generate_minimap2_pre_combination_event_dict " + 
                         "of multi_species_events.py.  Exiting . . . ")

    gtf_file.close()

    return pre_combination_event_dict


def evaluate_minimap_event_combinations(pre_combination_event_dict):
    '''
        Evaluates possible combinations of minimapped included and
        excluded event isoforms (i.e. if there is more than one of
        either) and performs preliminary check on whether they
        make sense.
    '''

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    used_ids = {}

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

    for event, event_val in pre_combination_event_dict.iteritems():

        event_type = event_val["event_type"]

        if len(event_val["included_forms"]) == 0:

            if len(event_val["excluded_forms"]) == 0: ######Is this a possible condition at all?  If not, remove it!!!

                continue

            else:

                if event_type in ["MR","RI"]:

                    for exc_version in event_val["excluded_forms"]:

                        exc_chrom = event_val["excluded_forms"][exc_version]["chrom"]
                        exc_strand = event_val["excluded_forms"][exc_version]["strand"]


                        exc_exons = event_val["excluded_forms"][exc_version]["exons"]

                        if ((event_type == "RI" and len(exc_exons) == 2) or 
                            (event_type == "MR" and len(exc_exons) >= 3)):

                            inc_exons = [[event_val["excluded_forms"][exc_version]["exons"][0][0], 
                                          event_val["excluded_forms"][exc_version]["exons"][1][1]]]

                            ######This is actually eliminating possibilities.  By using "event" as the indexing variable, any subsequent excluded form versions will overwrite this one.

                            if event in used_ids:

                                used_ids[event] += 1

                                new_id = (event + 
                                          "&" + 
                                          alphabet[used_ids[event]])

                            else:

                                used_ids[event] = 0
                                new_id = (event + 
                                          "&" + 
                                          alphabet[used_ids[event]])                                

                            minimap2_standard_event_dict.setdefault(event, {}
                                ).setdefault(new_id, {
                                    "event_type": event_type,
                                    "included_exons": copy.deepcopy(inc_exons),
                                    "excluded_exons": copy.deepcopy(exc_exons),
                                    "chrom": exc_chrom,
                                    "strand": exc_strand
                                })


        elif len(event_val["excluded_forms"]) == 0:

            if event_type in ["SE", "MS"]:


                for inc_version in event_val["included_forms"]:

                    inc_chrom = event_val["included_forms"][inc_version]["chrom"]
                    inc_strand = event_val["included_forms"][inc_version]["strand"]

                    inc_exons = event_val["included_forms"][inc_version]["exons"]

                    if (event_type is "SE" and len(inc_exons) == 3) or (event_type is "MS" and len(inc_exons) > 3):

                        exc_exons = [[event_val["included_forms"][inc_version]["exons"][0][0], 
                                      event_val["included_forms"][inc_version]["exons"][0][1]], 
                                     [event_val["included_forms"][inc_version]["exons"][-1][0], 
                                      event_val["included_forms"][inc_version]["exons"][-1][1]]]

                        ######Same situation as noted on line 183 - using event as index will cause subsequent included forms to overwrite this one.

                        if event in used_ids:

                            used_ids[event] += 1

                            new_id = (event + 
                                      "&" + 
                                      alphabet[used_ids[event]])

                        else:

                            used_ids[event] = 0
                            new_id = (event + 
                                      "&" + 
                                      alphabet[used_ids[event]])                                

                        minimap2_standard_event_dict.setdefault(event, {}
                            ).setdefault(new_id, {
                                "event_type": event_type,
                                "included_exons": copy.deepcopy(inc_exons),
                                "excluded_exons": copy.deepcopy(exc_exons),
                                "chrom": exc_chrom,
                                "strand": exc_strand
                            })                    



        else:

            for inc_version in event_val["included_forms"]:

                inc_chrom = event_val["included_forms"][inc_version]["chrom"]
                inc_strand = event_val["included_forms"][inc_version]["strand"]

                inc_exons = event_val["included_forms"][inc_version]["exons"]

                for exc_version in event_val["excluded_forms"]:

                    exc_chrom = event_val["excluded_forms"][exc_version]["chrom"]
                    exc_strand = event_val["excluded_forms"][exc_version]["strand"]

                    exc_exons = event_val["excluded_forms"][exc_version]["exons"]

                    def add_to_dict():

                        if event in used_ids:

                            used_ids[event] += 1

                            new_id = (event + 
                                      "&" + 
                                      alphabet[used_ids[event]])

                        else:

                            used_ids[event] = 0
                            new_id = (event + 
                                      "&" + 
                                      alphabet[used_ids[event]])                                

                        minimap2_standard_event_dict.setdefault(event, {}
                            ).setdefault(new_id, {
                                "event_type": event_type,
                                "included_exons": copy.deepcopy(inc_exons),
                                "excluded_exons": copy.deepcopy(exc_exons),
                                "chrom": exc_chrom,
                                "strand": exc_strand
                            })   

                    if inc_chrom == exc_chrom and inc_strand == exc_strand:

                        if event_type == "MS":

                            if (len(inc_exons) >= event_type_length_dict[event_type]["inc"] and 
                                len(exc_exons) == event_type_length_dict[event_type]["exc"]):

                                add_to_dict()

                        elif event_type == "MR":

                            if (len(inc_exons) == event_type_length_dict[event_type]["inc"] and 
                                len(exc_exons) >= event_type_length_dict[event_type]["exc"]):

                                add_to_dict()

                        elif event_type in ["SE", "RI", "A3", "A5", "AF", "AL", "MX"]:

                            if (len(inc_exons) == event_type_length_dict[event_type]["inc"] and 
                                len(exc_exons) == event_type_length_dict[event_type]["exc"]):

                                add_to_dict()

                        elif event_type in ["CO", "CF", "CL"]:

                            if len(inc_exons) >= 2 and len(exc_exons) >= 2:

                                add_to_dict()

                        elif event_type in ["MF", "ML"]:

                            if ((len(inc_exons) >= 3 and len(exc_exons) >= 2) or 
                                (len(inc_exons) >= 2 and len(exc_exons) >= 3)):

                                add_to_dict()


    return minimap2_standard_event_dict



def resolve_minimap2_event_combinations(
        minimap2_standard_event_dict):

    final_dict = {}

    mod_id_original_dict = {}

    for event, event_val in minimap2_standard_event_dict.iteritems():

        if len(event_val) == 1:

            final_dict[event] = event_val.values()[0]

        else:

            for new_id, new_id_val in event_val.iteritems():

                final_dict[new_id] = new_id_val

                mod_id_original_dict[new_id] = event

    return final_dict, mod_id_original_dict



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

                    new_ID = (ID.split("_")[0] + 
                              "&" + 
                              alphabet[used_IDs[ID]] + 
                              "_" + 
                              ID.split("_")[1])

                    entry[0] = new_ID

                else:

                    used_IDs[ID] = 0

                    new_ID = (ID.split("_")[0] + 
                              "&" + 
                              alphabet[used_IDs[ID]] + 
                              "_" + 
                              ID.split("_")[1])

                    entry[0] = new_ID

                new_samfile.write("\t".join(entry) + "\n")

            else:

                new_samfile.write(line.strip() + "\n")

    new_samfile.close()



def get_cross_species_events(
        from_species, 
        to_species, 
        species_dict, 
        num_threads, 
        minimap2_path = "minimap2", 
        samtools_path = "samtools", 
        bedtools_path = "bedtools", 
        bedToGenePred_path = "bedToGenePred", 
        genePredToGtf_path = "genePredToGtf"):

    #/hive/users/anjowall/minimap2/minimap2 -a -G 500000 -I 100G -t 10 -u f -x splice /hive/users/anjowall/genomes/rheMac8/STAR_genome_rheMac8/rheMac8.fa splice_lib_events.fa > minimap2_test/hu_to_rh_test.sam

    sam_outfile_path = (species_dict[to_species]["minimap_outdir"] + 
                        "/" + 
                        to_species + 
                        "_" + 
                        "minimapped_from_" + 
                        from_species + 
                        "_splice_events.sam")

    newsam_path = (sam_outfile_path.split(".sam")[0] + 
                   ".dedup_id.sam")


    try: ##from Paul Finn http://www.pfinn.net/python-check-if-file-exists.html
        with open(newsam_path[:-3] + "gtf") as file:

            print newsam_path[:-3] + "gtf exists - skipping minimap run."

    except IOError:

        try:

            subprocess.check_output(minimap2_path + 
                                   " -a -G 300000 --secondary=no -t " + 
                                   num_threads + 
                                   " -u f -x splice " + 
                                   species_dict[to_species]["genome_fasta"] + 
                                   " " + 
                                   species_dict[from_species]["event_fasta_path"] + 
                                   " > " + 
                                   sam_outfile_path, 
                                shell = True, 
                                stderr = subprocess.STDOUT)

            deduplicated_minimap2_sam(sam_outfile_path, newsam_path)

            subprocess.check_output(samtools_path + 
                                    " view -hb " + 
                                    newsam_path + 
                                    " > " + 
                                    newsam_path[:-3] + 
                                    "bam", 
                                shell = True, 
                                stderr = subprocess.STDOUT)

            subprocess.check_output(samtools_path + 
                                    " sort " + 
                                    newsam_path[:-3] + 
                                    "bam" + 
                                    " > " + 
                                    newsam_path[:-3] + 
                                    "sorted.bam", 
                                shell = True, 
                                stderr = subprocess.STDOUT)

            subprocess.check_output(samtools_path + 
                                    " index " + 
                                    newsam_path[:-3] + 
                                    "sorted.bam", 
                                shell = True, 
                                stderr = subprocess.STDOUT)

            subprocess.check_output(bedtools_path + 
                                    " bamtobed -bed12 -i " + 
                                    newsam_path[:-3] + 
                                    "sorted.bam" + 
                                    " > " + 
                                    newsam_path[:-3] + 
                                    "bed12", 
                                shell = True, 
                                stderr = subprocess.STDOUT)

            subprocess.check_output(bedToGenePred_path + 
                                    " " + 
                                    newsam_path[:-3] + 
                                    "bed12 stdout | " + 
                                    genePredToGtf_path + 
                                    " file stdin stdout | awk '$3==\"exon\"' > " + 
                                    newsam_path[:-3] + 
                                    "gtf", 
                                shell = True, 
                                stderr = subprocess.STDOUT)

        ##import as splice_lib event dict --- this needs something special to parse modified event id


        ###Need some code to deal with orphaned forms - e.g. orphaned SE included form actually serves as entire form etc etc (evaluate_orphans function)

        ###Also run some complete event dict code

        except subprocess.CalledProcessError as e:

            print e.output

            sys.exit("minimap2 failed in " + 
                     "'get_cross_species_de_novo_events' " + 
                     "of multi_species_events_standalone.py. " + 
                     "Exiting . . . ")

    precombination_event_dict = generate_minimap2_pre_combination_event_dict(
        newsam_path[:-3] + "gtf")
    minimap2_standard_event_dict = evaluate_minimap_event_combinations(
        precombination_event_dict)

    minimap2_standard_event_dict, id_translator = resolve_minimap2_event_combinations(
        minimap2_standard_event_dict)

    return minimap2_standard_event_dict, id_translator





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


    if ((initial_event_type in ["AF","MF","CF","UF"] and strand == "+") or 
        (initial_event_type in ["AL","ML","CL","UL"] and strand == "-")):

        if ((event_dict_entry["included_exons"][0][0] != 
             event_dict_entry["excluded_exons"][0][0]) and 
            (event_dict_entry["included_exons"][-1][-1] == 
             event_dict_entry["excluded_exons"][-1][-1])):

            form_1_exons.insert(0, "universal_left")
            form_2_exons.insert(0, "universal_left")
            pre_common.append("universal_left")

        else:

            #print event, "Pre-comprehensive misclassification tag"
            event_misclassified = True ### replace this with something more helpful


    elif ((initial_event_type in ["AF","MF","CF","UF"] and strand == "-") or 
          (initial_event_type in ["AL","ML","CL","UL"] and strand == "+")):

        if ((event_dict_entry["included_exons"][-1][-1] != 
             event_dict_entry["excluded_exons"][-1][-1]) and 
            (event_dict_entry["included_exons"][0][0] == 
             event_dict_entry["excluded_exons"][0][0])):

            form_1_exons.append("universal_right")
            form_2_exons.append("universal_right")
            post_common.append("universal_right")


        else:

            event_misclassified = True ### replace this with something more helpful
            #print event, "Pre-comprehensive misclassification tag"

    else:

        if not ((event_dict_entry["included_exons"][0][0] == 
                 event_dict_entry["excluded_exons"][0][0]) and 
                (event_dict_entry["included_exons"][-1][-1] == 
                 event_dict_entry["excluded_exons"][-1][-1])):

            event_misclassified = True ### replace this with something more helpful
            #print event, "Pre-comprehensive misclassification tag"

    if not event_misclassified:

        class_result = infer_pairwise_events.classify_event(
            pre_common, 
            post_common, 
            tx1_unique, 
            tx2_unique, 
            form_1_exons, 
            form_2_exons, 
            strand, 
            tx1, 
            tx2)

        #print event, class_result

        same_classification = initial_event_type == class_result["event_type"]

        return same_classification

    else:

        return False


def get_junction_key(event_dict_entry):

    flat_joined_included_jl = "_".join([str(i) for j in 
                                        event_dict_entry["included_junctions"] 
                                        for i in j])

    flat_joined_excluded_jl = "_".join([str(i) for j in 
                                        event_dict_entry["excluded_junctions"] 
                                        for i in j])

    flat_joined_included_eij = "_".join(map(
        str, 
        event_dict_entry["included_ei_junctions"]))
    flat_joined_excluded_eij = "_".join(map(
        str, 
        event_dict_entry["excluded_ei_junctions"]))


    chrom = event_dict_entry["chrom"]
    strand = event_dict_entry["strand"]

    key = (chrom + 
           "_" + 
           flat_joined_included_jl + 
           "|" + 
           flat_joined_included_eij + 
           ":" +  
           flat_joined_excluded_jl + 
           "|" + 
           flat_joined_excluded_eij + 
           "_" + 
           strand)


    return key


def merge_event_dicts_by_species(
        event_dicts, 
        cross_species_event_dicts, 
        species_in_order, 
        outdir):

    merged_event_dict = add_native_events_to_merged_dict(
                            event_dicts)

    add_cross_species_events_to_merged_dict(
        merged_event_dict, 
        event_dicts,
        cross_species_event_dicts)

    del event_dicts
    del cross_species_event_dicts

    reassign_to_native_only_list = reconcile_cross_species_information(
        merged_event_dict,
        species_in_order)

    remove_ambiguously_mapped_events(
        reassign_to_native_only_list, 
        merged_event_dict)

    return merged_event_dict



def add_native_events_to_merged_dict(event_dicts):

    merged_events = {}

    for species, species_val in event_dicts.iteritems():

        print "event_dict for species", species, "has length", len(species_val)

        merged_events[species] = {}

        for event, event_val in species_val.iteritems():

            junction_key = get_junction_key(event_val)

            merged_events[species][junction_key] = copy.deepcopy(event_val)
            new_entry = merged_events[species][junction_key]

            new_entry["native"] = True
            new_entry["native_id"] = event
            new_entry["cross_species_keys"] = {}
            new_entry["cross_species_keys"][species] = {event: junction_key}

    return merged_events



def add_cross_species_events_to_merged_dict(
        merged_events, 
        event_dicts,
        cross_species_event_dicts):

    for species, species_val in cross_species_event_dicts.iteritems():

        for other_species, other_species_val in species_val.iteritems():

            print ("cross_species_event_dict for species " +
                   species + 
                   " from other_species " + 
                   other_species + 
                   " has length " +
                   str(len(other_species_val)))

            for event, event_val in other_species_val.iteritems():


                junction_key = get_junction_key(
                        event_val)


                if event in event_dicts[other_species]:

                    other_junction_key = get_junction_key(
                        event_dicts[other_species][event])  
                        ## get the "native" event junction key

                else:

                    print >> sys.stderr, event, "Not found in event dict of", other_species
                    other_junction_key = None

                if junction_key in merged_events[species]:

                    merged_event_entry = merged_events[species][junction_key]["cross_species_keys"]

                    species_event_id = merged_event_entry[species].keys()[0]
                    ###Must reference the junction key from the source species !!!!!!

                    if other_junction_key:
                        # the native event will have already been added to the merged dict

                        merged_event_entry[other_species] = {event: other_junction_key}

                        merged_event_entry_other = merged_events[other_species][other_junction_key]["cross_species_keys"]

                        merged_event_entry_other[species] = {species_event_id: junction_key}

                else:

                    merged_events[species][junction_key] = copy.deepcopy(event_val)
                    merged_event_new_entry = merged_events[species][junction_key]
                    merged_event_new_entry["native"] = False
                    merged_event_new_entry["native_id"] = None

                    if other_junction_key:

                        merged_event_new_entry["cross_species_keys"] = {
                                        other_species: {
                                                event: other_junction_key
                                                }, 
                                              species: {
                                                "placeholder": junction_key
                                                }
                                            }

                        merged_event_entry_other = merged_events[other_species][other_junction_key]["cross_species_keys"]

                        merged_event_entry_other[species] = {"placeholder": junction_key}

                    else:

                        merged_event_new_entry["cross_species_keys"] = {species: {"placeholder": junction_key}}




def reconcile_cross_species_information(
        merged_events,
        species_in_order):

    ## At this point all native events should have all of the information about where they successfully cross-mapped

    reassign_to_native_only_list = []

    for species, species_val in merged_events.iteritems():

        for event, event_entry in species_val.iteritems():

            #print "event"
            #print event

            master_cross_species_keys = copy.deepcopy(event_entry["cross_species_keys"])

            ## Now go back through the cross-species keys and extract the cross-species key information coming from the other species dicts
            ## If repeated information (i.e. about a specific species) is found, check that it matches what is already present.  If it is not already present, add it.

            reassign_to_native_only = False

            for cross_species, cross_species_val in event_entry["cross_species_keys"].iteritems():

                for cross_id, cross_junction in cross_species_val.iteritems():

                    try:

                        for check_species, check_species_entry in merged_events[cross_species][cross_junction]["cross_species_keys"].iteritems():

                            if check_species in master_cross_species_keys:

                                if check_species_entry != master_cross_species_keys[check_species]:

                                    # Do the following instead:
                                    # 1) break this relationship
                                    # 2) if events in this comparison are solely defined by                          

                                    reassign_to_native_only = True
                                    reassign_to_native_only_list.append(master_cross_species_keys)
                                    reassign_to_native_only_list.append(
                                        copy.deepcopy(merged_events[cross_species][cross_junction]["cross_species_keys"]))


                            else:

                                master_cross_species_keys.setdefault(
                                    check_species, 
                                    copy.deepcopy(check_species_entry))

                    except KeyError:

                        print " "
                        print cross_species
                        print "_"
                        print cross_junction
                        print " "
                        print event_entry
                        pickle.dump(merged_events, open("./merged_events_pre_rename_during_rec_standard_event_dict.pkl", 'wb'), pickle.HIGHEST_PROTOCOL)
                        sys.exit()

            if not reassign_to_native_only:

                for reassign_species, reassign_entry in master_cross_species_keys.iteritems():

                    reassign_key = reassign_entry.values()[0]

                    merged_events[reassign_species][reassign_key]["cross_species_keys"] = copy.deepcopy(master_cross_species_keys)


    return reassign_to_native_only_list



def remove_ambiguously_mapped_events(
        reassign_to_native_only_list, 
        merged_events):

    

    non_native_keys_to_delete = {}

    for key_dict in reassign_to_native_only_list:

        for species, species_entry in key_dict.iteritems():

            key = species_entry.values()[0]


            merged_events_entry = merged_events[species][key]

            if merged_events_entry["native"]:

                merged_events_entry["ambiguous_mapping_issue"] = True

                for species_key in merged_events_entry["cross_species_keys"].keys():

                    if species_key != species:

                        del merged_events_entry["cross_species_keys"][species_key]

            else:
                non_native_keys_to_delete.setdefault(species, set()).add(key)

    ## delete non-native ambiguously mapping keys

    for species, key_set in non_native_keys_to_delete.iteritems():

        for key in key_set:

            del merged_events[species][key]


def rename_merged_events(
        merged_events, 
        species_in_order):

    #pickle.dump(merged_events, open(outdir + "/merged_events_pre_rename_standard_event_dict.pkl", 'wb'), pickle.HIGHEST_PROTOCOL)

    #human
    #macaque
    #chr7_|49832632_49831973:49831972_49832633|_-
    #{'excluded_ei_junctions': [], 'event_type': 'RI', 'sources': [], 'native_id': 'RI.0001510', 'included_exons': [['74034478', '74036394']], 'included_junction_counts': {'chr15_74035171_+': 0, 'chr15_74034531_+': 0}, 'included_form_transcripts': [], 'included_count': 0, 'excluded_exons': [['74034478', '74034530'], ['74035172', '74036394']], 'native': True, 'excluded_jn_count': '3', 'included_junctions': [], 'included_unique_edges': [], 'excluded_count': 0, 'excluded_junction_counts': {'chr15_left_74035172_+': 0, 'chr15_right_74034530_+': 0, 'chr15_74034530_74035172_+': 0}, 'excluded_junctions': [['74034530', '74035172']], 'included_jn_count': '2', 'chrom': 'chr15', 'strand': '+', 'included_ei_junctions': [74034531, 74035171], 'excluded_unique_edges': ['left_74035172', 'right_74034530'], 'excluded_form_transcripts': [], 'cross_species_keys': {'macaque': {'RI.0002500': 'chr7_|49832632_49831973:49831972_49832633|_-'}}, 'gene': ['PML']}
    #chr15_|74034531_74035171:74034530_74035172|_+



    renamed_merged_events = {}

    event_counter = 1

    for species in species_in_order:

        merged_species = merged_events[species]

        while len(merged_events[species]) > 0:

            working_event_key = merged_species.keys()[0]
            #print working_event_key
            working_event = merged_species[working_event_key]
            working_event_cs_keys = working_event["cross_species_keys"]

            if species not in working_event_cs_keys:

                sys.exit("Native species not in cs keys")

            species_in_id = [ref_species for ref_species 
                             in species_in_order 
                             if ref_species in 
                             working_event_cs_keys]

            all_species_keys = {ref_species: working_event_cs_keys[ref_species].values()[0] 
                                for ref_species in species_in_id}


            new_id_primary = (working_event["event_type"] + 
                              "." + 
                              str(event_counter).zfill(7))

            full_new_id = ("_".join(species_in_id) + 
                           "-" + 
                           new_id_primary)

            for ref_species in all_species_keys:

                ref_species_event_entry = (merged_events[ref_species]
                                            ).get(all_species_keys[ref_species])

                if ref_species_event_entry:

                    ((renamed_merged_events
                        ).setdefault(ref_species, {})
                            ).setdefault(
                                full_new_id, 
                                ref_species_event_entry)

                    del merged_events[ref_species][all_species_keys[ref_species]]

            event_counter += 1

    return renamed_merged_events



def complete_merged_dict(
        merged_event_dict, 
        species_dict, 
        add_gene_symbols, 
        transcript_gtfs):

    ##add junctions, junction counts slots, remove duplicates or otherwise problematic events, clean up outermost exon edges


    for species, species_val in merged_event_dict.iteritems():

        splice_lib.complete_event_dict(
            species_val,
            inform_using_ri_events = False,
            suppress_unique_edges = True)

        species_dict_entry = species_dict[species]
        annotated_transcript_gtf = species_dict_entry["annotated_transcript_gtf"]

        if add_gene_symbols:

            if transcript_gtfs:

                (matches_dict, 
                 genes_dict, 
                 transcript_dict) = assign_events_to_genes.main(
                        ["--transcript_gtf", 
                         species_dict_entry["merged_stringtie_gtf_path"], 
                         "--gene_gtf", 
                         annotated_transcript_gtf, 
                         "--prefix", 
                         "placeholder", 
                         "--suppress_output"], 
                         species_val)

            else:

                (matches_dict, 
                 genes_dict, 
                 transcript_dict) = assign_events_to_genes.main(
                        ["--gene_gtf", 
                         annotated_transcript_gtf, 
                         "--prefix", 
                         "placeholder", 
                         "--suppress_output"], 
                         species_val)


            for gene, gene_val in matches_dict.iteritems():

                for event in gene_val:

                    (species_val[event]
                        ).setdefault("gene", []).append(gene)




def check_gene_equivalence(
        merged_event_dict, 
        outdir, 
        species_in_order):

    '''
        Check that an event maps to the equivalent gene (using CAT annotations or similar) in all species.  If it doesn't, write the information to an output file so it can be investigated and/or discarded in downstream analysis
    '''

    events_genes = {}

    multiple_genes = open(outdir + "/events_with_contrary_genes.tsv", 'w')

    multiple_genes.write("\t".join(["event_id"] + 
                         species_in_order + 
                         ["consensus_overlap"]) + 
                         "\n")


    for species, species_val in merged_event_dict.iteritems():

        for event, event_val in species_val.iteritems():

            if "gene" in event_val:

                ((events_genes
                    ).setdefault(event, {})
                        ).update({species: ",".join(event_val["gene"])})


    for event, event_val in events_genes.iteritems():

        all_genes = []

        for species in species_in_order:

            if species in event_val:

                all_genes.append(event_val[species])


        if len(set(all_genes)) > 1:

            consensus_overlap = False

            if len(set.intersection(*[set(i.split(",")) 
                                    for i in all_genes])) > 0:

                consensus_overlap = True

            genes_to_write = []

            for species in species_in_order:

                if species in event_val:

                    genes_to_write.append(event_val[species])

                else:

                    genes_to_write.append("NA")

            multiple_genes.write("\t".join([event] + 
                                 genes_to_write + 
                                 [str(consensus_overlap)]) + 
                                 "\n")

    multiple_genes.close()



def run_paftools_liftover(
        path_to_paf_file,
        input_gtf_path,
        outdir,
        from_species,
        to_species,
        paftools_path,
        min_length = 5000):

    gtf_dict = {}

    with open(input_gtf_path, 'r') as file:

        for line in file:

            line_strip = line.strip()

            entry = line.split("\t")

            key = ";".join( [ entry[0], str(int(entry[3]) - 1), str(entry[4]) ] )


            gtf_dict.setdefault(key, set()).add(line_strip)



    paftools_bed_out = outdir + "/" + from_species + "_" + to_species + "_paftools.bed"

    paftools_gtf_out = outdir + "/" + from_species + "_" + to_species + "_paftools.gtf"


    paftools_lift_cmd = (paftools_path + 
                         " liftover " + 
                         "-l " + 
                         str(min_length) + 
                         " " + 
                         path_to_paf_file + 
                         ' <(awk \'{OFS="\\t"; print $1,$4-1,$5}\' ' + 
                         input_gtf_path + " | sort | uniq ) > " + 
                         paftools_bed_out) 

    try: ##from Paul Finn http://www.pfinn.net/python-check-if-file-exists.html
        with open(paftools_bed_out) as file:
            print (paftools_bed_out + 
                   " bed exists - skipping paftools liftover run.")
    except IOError:
        try:

            print "Running"
            print paftools_lift_cmd
            gm.run_bash_cmd(paftools_lift_cmd.split())
        except subprocess.CalledProcessError as e:
            print e.output
            sys.exit("paftools liftover failed in " + 
                     "'run_paftools_liftover' of " + 
                     "multi_species_events_standalone.py. " + 
                     "Exiting . . . ")    

    old_new_interval_dict = {}

    with open(paftools_bed_out, 'r') as file:

        for line in file:

            entry = line.split()
            old_key = ";".join(entry[3].split(";")[0:3])
            new_value = ";".join( [ entry[0], entry[1], entry[2], entry[5] ] )

            old_new_interval_dict.setdefault(
                old_key, set()).add(new_value)


    with open(paftools_gtf_out, 'w') as file:

        for old, new in old_new_interval_dict.iteritems():

            if len(new) == 1:

                new_value = new.pop()

                if old in gtf_dict:

                    for gtf_entry in gtf_dict[old]:

                        new_gtf_entry = gtf_entry.split("\t")

                        new_split = new_value.split(";")

                        new_gtf_entry[0] = new_split[0]
                        new_gtf_entry[3] = str(int(new_split[1]) + 1)
                        new_gtf_entry[4] = new_split[2]

                        if new_split[3] == "-":

                            if new_gtf_entry[6] == "+":

                                new_gtf_entry[6] = "-"

                            elif new_gtf_entry[6] == "-":

                                new_gtf_entry[6] = "+"

                        file.write("\t".join(new_gtf_entry) + "\n")

    return paftools_gtf_out


def run_liftover(
        path_to_chain_file, 
        input_gtf_path, 
        outdir, 
        from_species, 
        to_species,
        liftover_path):

    path_to_successful_lifts = (outdir + 
                                "/" + 
                                from_species + 
                                "_to_" + 
                                to_species + 
                                "_lift.gtf")

    try: ##from Paul Finn http://www.pfinn.net/python-check-if-file-exists.html
        with open(path_to_successful_lifts) as file:

            print (path_to_successful_lifts + 
                   " gtf exists - skipping liftOver run.")

    except IOError:

        try:

            cmd_str = (liftover_path + 
                       " minMatch=0.8 -gff " + 
                       input_gtf_path + 
                       " " + 
                       path_to_chain_file + 
                       " " + 
                       path_to_successful_lifts + 
                       " " + 
                       outdir + 
                       "/" + 
                       from_species + 
                       "_to_" + 
                       to_species + 
                       "_lift_fail.gtf")

            gm.run_bash_cmd(cmd_str.split())


        except subprocess.CalledProcessError as e:

            print e.output

            sys.exit("liftOver failed in " + 
                     "'run_liftover' of " + 
                     "multi_species_events_standalone.py. " + 
                     "Exiting . . . ")

    return path_to_successful_lifts



def same_exon_counts(event_dict_entry_list):
    '''
        takes as input a list of event dict entries 
        and checks whether they all have the same number 
        of included and excluded form exons.  Returns boolean.
    '''

    included_exon_counts = [len(i["included_exons"]) 
                            for i in 
                            event_dict_entry_list]

    excluded_exon_counts = [len(i["excluded_exons"]) 
                            for i in event_dict_entry_list]

    if len(set(included_exon_counts)) > 1 or len(set(excluded_exon_counts)) > 1:

        return False

    return True



def check_exon_intron_length_overlap(
        event_entry,
        min_exon_length = 3,
        min_intron_length = 20):


    violation_list = set()

    included_exons = event_entry["included_exons"]
    excluded_exons = event_entry["excluded_exons"]

    included_junctions = event_entry["included_junctions"]
    excluded_junctions = event_entry["excluded_junctions"]


    for exon in included_exons + excluded_exons:
            if (int(exon[1]) - int(exon[0]) + 1) < min_exon_length:
                violation_list.add("exon_length_violation")
                break

    for form in ["included", "excluded"]:
            for exon in event_entry[form + "_exons"]:
                for other_exon in event_entry[form + "_exons"]:
                    if exon != other_exon:
                        if (int(exon[0]) in range(int(other_exon[0]), int(other_exon[1]) + 1) or 
                            int(exon[1]) in range(int(other_exon[0]), int(other_exon[1]) + 1)):
                            violation_list.add("exon_overlap_violation")
                            break

    for junction in included_junctions + excluded_junctions:
            if (int(junction[1]) - int(junction[0]) + 3) < min_intron_length:
                violation_list.add("intron_length_violation")
                break        

    if len(included_exons) == 0 or len(excluded_exons) == 0:

        violation_list.add("missing_form_violation")

    else:

        if (not (int(included_exons[0][0]) >= int(excluded_exons[0][0]) and 
            int(included_exons[0][0]) <= int(excluded_exons[-1][-1]) + 1) and 
            not (int(included_exons[-1][-1]) >= int(excluded_exons[0][0]) and 
            int(included_exons[-1][-1]) <= int(excluded_exons[-1][-1]) + 1) and 
            not (int(excluded_exons[0][0]) >= int(included_exons[0][0]) and 
            int(excluded_exons[0][0]) <= int(included_exons[-1][-1]) + 1) and 
            not (int(excluded_exons[-1][-1]) >= int(included_exons[0][0]) and 
            int(excluded_exons[-1][-1]) <= int(included_exons[-1][-1]) + 1)):
            
            violation_list.add("form_overlap_violation")

    return violation_list











def check_mapped_event_integrity(
        mapped_species, 
        source_species, 
        mapped_events_dict, 
        source_events_dict, 
        mapping_failure_log_file,
        id_translator = None):
    '''
        Ensures that mapped events have same exon counts for both source and mapped features, as well as the same event type

        Assumes that mapped_events_dict and source_events_dict have same keys
    '''


    event_issues = {}

    for event, event_val in mapped_events_dict.iteritems():

        write_to_log = False

        violations = check_exon_intron_length_overlap(
            event_val)

        if len(violations) > 0:

            _ = [event_issues.setdefault(event, set()).add(i) for i in violations]
            write_to_log = True

        try:

            if not same_exon_counts([event_val, 
                                     source_events_dict[event]]):

                event_issues.setdefault(event, set()).add("exon_count")
                write_to_log = True

        except KeyError:

            if id_translator:

                try:

                    if not same_exon_counts([event_val, 
                                             source_events_dict[id_translator[event]]]):

                        event_issues.setdefault(event, set()).add("exon_count")
                        write_to_log = True

                except KeyError:

                    print event
                    print event_val
                    sys.exit(
                        "Problem finding event", 
                        event, 
                        " in source dict in check_mapped_event_integrity()")

            else:

                sys.exit("Problem finding event", 
                         event, 
                         " in source dict in check_mapped_event_integrity()")



        if not (len(event_val["included_exons"]) == 0 or 
                len(event_val["excluded_exons"]) == 0):

            #print mapped_events_dict[event]["included_exons"]
            #print mapped_events_dict[event]["excluded_exons"]

            same_classification = check_event_type(event_val, event)

            if not same_classification:

                event_issues.setdefault(event, set()
                    ).add("post_map_class_error")

                write_to_log = True

        else:

            same_classification = False

            event_issues.setdefault(event, set()
                ).add("post_map_class_error")

            write_to_log = True

        if write_to_log:

            mapping_failure_log_file.write(
                "\t".join([mapped_species + 
                           "," + 
                           source_species, 
                           event, 
                           ",".join(event_issues[event])]) + 
                "\n")


    ### Now actually remove the offending event


    if id_translator:

        for event in event_issues:

            del mapped_events_dict[event]

            if event in id_translator:

                del id_translator[event]

        return id_translator

    else:

        for event in event_issues:

            del mapped_events_dict[event]



def assess_multiply_mapped_events(
        id_translator,
        minimap2_standard_event_dict,
        source_species,
        target_species,
        multimapper_log):
    '''
    Assess whether any multiply mapped events remain so
    after pruning based on integrity checks (i.e. some
    possibilities may be removed due to inconsistencies
    with original event).  Removes persistently multimapped
    events from event dict and writes them to log file.

    Parameters
    ----------

    id_translator : dict
        dict with new ids (i.e. original event ids with 
        alphabetic suffix to discriminate amongst multiple
        possibilities) as keys and original ids as values
    minimap2_standard_event_dict : dict
        event dict derived from reading-in and processing
        minimap2 event mappings
    source_species : str
        Source species from which the mapped event sequences
        were originally derived
    target_species : str
        Target species to which the original mapped event
        sequences were mapped
    multimapper_log : file
        Log file containing information about persistently
        multimapped events.
    '''

    inverted_translator = {}

    for new_id, event in id_translator.iteritems():

        inverted_translator.setdefault(
            event, set()).add(new_id)

    for event, new_ids in inverted_translator.iteritems():

        if len(new_ids) == 1:

            new_id = list(new_ids)[0]

            new_id_val = copy.deepcopy(
                minimap2_standard_event_dict[new_id])

            del minimap2_standard_event_dict[new_id]

            minimap2_standard_event_dict[event] = new_id_val

        else:

            multimapper_log_entry = "\t".join(
                [source_species, 
                 target_species,
                 event]) + "\n"

            multimapper_log.write(multimapper_log_entry)



def final_original_name(
    merged_event_dict, 
    species_list, 
    outdir):

    #merged_events[species][junction_key]["cross_species_keys"][species]


    species_comparisons = []

    for i, species_1 in enumerate(species_list[:-1]):

        for j, species_2 in enumerate(species_list[i+1:]):

            species_comparisons.append([species_1,species_2])

    entry_set = set()

    for species, species_val in merged_event_dict.iteritems():

        for event, event_val in species_val.iteritems():

            event_val_csk = event_val["cross_species_keys"]

            for comparison in species_comparisons:

                comparison_0 = comparison[0]
                comparison_1 = comparison[1]

                if (comparison_0 in event_val_csk and 
                    comparison_1 in event_val_csk):

                    event_val_csk_0 = event_val_csk[comparison_0].keys()
                    event_val_csk_1 = event_val_csk[comparison_1].keys()


                    if (len(event_val_csk_0) > 1 or 
                       len(event_val_csk_1) > 1):

                        print ("WARNING: event " +
                               event + 
                               " has more than one " + 
                               "cross-species key for " +
                               "species " + comparison_0 +
                               " or " + comparison_1)

                    entry_string = "\t".join([
                        event, 
                        event_val_csk_0[0], 
                        event_val_csk_1[0], 
                        comparison_0, 
                        comparison_1])

                    entry_set.add(entry_string)

                elif comparison_0 in event_val_csk:

                    event_val_csk_0 = event_val_csk[comparison_0].keys()

                    if len(event_val_csk_0) > 1:

                        print ("WARNING: event " +
                               event + 
                               " has more than one " + 
                               "cross-species key for " +
                               "species " + comparison_0)


                    entry_string = "\t".join([
                        event, 
                        event_val_csk_0[0], 
                        "NA", 
                        comparison_0, 
                        "NA"])

                    entry_set.add(entry_string)

                elif comparison_1 in event_val_csk:

                    event_val_csk_1 = event_val_csk[comparison_1].keys()

                    if len(event_val_csk_1) > 1:

                        print ("WARNING: event " +
                               event + 
                               " has more than one " + 
                               "cross-species key for " +
                               "species " + comparison_1)


                    entry_string = "\t".join([
                        event, 
                        "NA", 
                        event_val_csk_1[0], 
                        "NA", 
                        comparison_1]) 

                    entry_set.add(entry_string)               


    with open(outdir + "/final_original_name.tsv", 'w') as file:

        file.write(
            "\t".join([
                "final_event_id", 
                "original_event_id_species_1", 
                "original_event_id_species_2", 
                "species_1", 
                "species_2"]) + "\n")

        for entry in entry_set:

            file.write(entry + "\n")




def write_output(
        renamed_merged_event_dict,
        outdir,
        species_list):

    check_gene_equivalence(
        renamed_merged_event_dict, 
        outdir, 
        species_list)


    for species, species_val in renamed_merged_event_dict.iteritems():

        splice_lib.output_event_gtf(
            species_val, 
            outdir, 
            name = species + "_splice_lib_events")

        splice_lib.output_event_bedfile(
            species_val, 
            outdir, 
            name = species + "_splice_lib_events")

        splice_lib.output_event_gtf_with_transcripts(
            species_val, 
            outdir, 
            name = species + "_splice_lib_events_with_transcripts")

        splice_lib.output_miso_event_gff3(
            species_val, 
            outdir, 
            name = species + "_splice_lib_events")


    final_original_name(
        renamed_merged_event_dict, 
        species_list, 
        outdir)



def open_logs(
        outdir,
        mode):

    mapping_failure_log = open(outdir + 
                               "/mapping_failure_log.tsv", 
                               'w')

    mapping_failure_log.write("\t".join(["mapped_species,source_species", 
                                         "event_id", 
                                         "event_issues"]) + 
                              "\n")

    if mode == "minimap":

        multimapper_log = open(outdir + "/multimapper_log.tsv", 'w')
        multimapper_log.write("\t".join([
            "source_species",
            "target_species",
            "event_id"]) + 
            "\n")

    else:

        multimapper_log = None

    return mapping_failure_log, multimapper_log



def ambiguous_mapper_log_write(
        renamed_merged_event_dict,
        outdir):

    with open(outdir + "/ambiguous_mapper_log.tsv", 'w') as file:

        file.write("\t".join([
            "source_species",
            "event_id"]) + 
            "\n")     

        native_events = {}

        for species, species_val in renamed_merged_event_dict.iteritems():

            for event, event_val in species_val.iteritems():

                if event_val["native"]:

                    native_events.setdefault(event, []
                        ).append(species)

                if event_val.get("ambiguous_mapping_issue"):

                    file.write("\t".join([species, event]) + "\n")

    return native_events



def write_native_events(
        native_events,
        outdir):

    with open(outdir + 
              "/native_event_species.tsv", "w") as file:

        file.write("\t".join(["event_id", 
                              "species_native"]) + 
                   "\n")

        for event, species_list in native_events.iteritems():

            file.write("\t".join([event, 
                                  ",".join(sorted(species_list))]) + 
                       "\n")


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

                match = transcript

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

    for event in final_matches.keys():

        if final_matches[event]["support"] != "both":

            del final_matches[event]

    return final_matches



def exon_count_cross_species_step(
        cross_species_event_dicts,
        other_species,
        species,
        ioe_paths_dict,
        species_dict,
        event_dicts):


    splice_lib.add_transcripts_to_event_dict(
        ioe_paths_dict[species], 
        event_dicts[species])

    splice_lib.add_transcripts_to_event_dict(
        ioe_paths_dict[other_species], 
        event_dicts[other_species])


    species_transcript_centric_included_forms = make_transcript_event_dict(
                                                        event_dicts[species], 
                                                        "included")
    species_transcript_centric_excluded_forms = make_transcript_event_dict(
                                                        event_dicts[species], 
                                                        "excluded")
    other_species_transcript_centric_included_forms = make_transcript_event_dict(
                                                        event_dicts[other_species], 
                                                        "included")
    other_species_transcript_centric_excluded_forms = make_transcript_event_dict(
                                                        event_dicts[other_species], 
                                                        "excluded")

    matches = find_matches(
        event_dicts[species], 
        event_dicts[other_species], 
        species_dict[species]["annotated_transcript_dict"], 
        species_dict[other_species]["annotated_transcript_dict"], 
        species_transcript_centric_included_forms, 
        species_transcript_centric_excluded_forms, 
        other_species_transcript_centric_included_forms, 
        other_species_transcript_centric_excluded_forms)


    for event, match in matches.iteritems():

         cross_species_event_dicts.setdefault(
            other_species, {}).setdefault(
                species, {}).setdefault(
                    event, copy.deepcopy(
                        event_dicts[other_species][match["match"]]))




def generate_ioe_dict(
        ioe_paths_dict,
        species):

    ioe_dict = {}

    with open(ioe_paths_dict[species], 'r') as file:

        next(file) # skip header

        for line in file:

            entry = line.split()

            included_transcripts = set(entry[3].split(","))

            all_transcripts = set(entry[4].split(","))

            excluded_transcripts = all_transcripts - included_transcripts

            for inc_tx in included_transcripts:

                for exc_tx in excluded_transcripts:

                    key = inc_tx + ":" + exc_tx

                    event_id = entry[2].split(";")[1]

                    ioe_dict.setdefault(key, set()).add(event_id)

    return ioe_dict



def ioe_files_cross_species_step(
        cross_species_event_dicts,
        other_species,
        species,
        ioe_paths_dict,
        event_dicts):

    other_species_ioe_dict = generate_ioe_dict(
        ioe_paths_dict,
        other_species)

    species_ioe_dict = generate_ioe_dict(
        ioe_paths_dict,
        species)

    for key, value in species_ioe_dict.iteritems():

        if (len(value) == 1 or 
            len({i.split(".")[0] for i in value}) == len(value)):
            ## the above ensures that either a) this transcript pair
            ## has only one differentiating AS event or b)
            ## the multiple differentiating AS events are of unique
            ## types

            for event in value:

                event_type = event.split(".")[0]
                other_value = other_species_ioe_dict.get(key)

                if other_value:

                    if (len(other_value) == 1 or 
                        len({i.split(".")[0] for i in other_value}) == len(other_value)):

                        for other_event in other_value:

                            other_event_type = other_event.split(".")[0]

                            if event_type == other_event_type:

                                cross_species_event_dicts.setdefault(
                                    other_species, {}).setdefault(
                                        species, {}).setdefault(
                                            event, copy.deepcopy(
                                                event_dicts[other_species][other_event]))






def paftools_lift_cross_species_step(
        liftover_cross_species_event_dicts,
        other_species,
        species,
        paf_paths_dict,
        species_dict,
        paftools_path,
        outdir):

    liftover_cross_species_event_dicts.setdefault(other_species, {}
        ).update({species: splice_lib.generate_standard_event_dict(
            run_paftools_liftover(
                        paf_paths_dict[species][other_species],
                        species_dict[species]["native_event_gtf_path"],
                        outdir,
                        species,
                        other_species,
                        paftools_path,
                        min_length = 5000))})


def liftover_cross_species_step(
        liftover_cross_species_event_dicts,
        other_species,
        species,
        chain_paths_dict,
        species_dict,
        liftover_path,
        outdir):

    liftover_cross_species_event_dicts.setdefault(other_species, {}
        ).update({species: splice_lib.generate_standard_event_dict(
            run_liftover(chain_paths_dict[species][other_species], 
                         species_dict[species]["native_event_gtf_path"], 
                         outdir, 
                         species, 
                         other_species,
                         liftover_path))})


def minimap_cross_species_step(
        species,
        other_species,
        species_dict,
        num_threads,
        minimap2_path,
        samtools_path,
        bedtools_path,
        bedToGenePred_path,
        genePredToGtf_path,
        cross_species_event_dicts):

    cross_dict, id_translator = get_cross_species_events(
            species, 
            other_species, 
            species_dict, 
            num_threads, 
            minimap2_path = minimap2_path, 
            samtools_path = samtools_path, 
            bedtools_path = bedtools_path, 
            bedToGenePred_path = bedToGenePred_path, 
            genePredToGtf_path = genePredToGtf_path)

    cross_species_event_dicts.setdefault(other_species, {}
        ).update({species: cross_dict})

    return cross_dict, id_translator






def run_cross_species_processes(
        species_dict,
        species_list,
        event_dicts,
        outdir,
        num_threads,
        minimap2_path,
        samtools_path,
        bedToGenePred_path,
        genePredToGtf_path,
        bedtools_path,
        add_gene_symbols,
        transcript_gtfs,
        mode,
        chain_paths_dict,
        liftover_path,
        paf_paths_dict,
        paftools_path,
        ioe_paths_dict):

    cross_species_event_dicts = {}

    mapping_failure_log, multimapper_log = open_logs(outdir, mode)

    for i,species in enumerate(species_list):

        for other_species in species_list:

            if other_species != species:

                id_translator = None

                if mode == "minimap":

                    print "Running minimap cross-species step " + species + " " + other_species

                    cross_dict, id_translator = minimap_cross_species_step(
                        species,
                        other_species,
                        species_dict,
                        num_threads,
                        minimap2_path,
                        samtools_path,
                        bedtools_path,
                        bedToGenePred_path,
                        genePredToGtf_path,
                        cross_species_event_dicts)

                elif mode == "liftover":

                    print "Running UCSC liftover cross-species step " + species + " " + other_species                    

                    liftover_cross_species_step(
                        cross_species_event_dicts,
                        other_species,
                        species,
                        chain_paths_dict,
                        species_dict,
                        liftover_path,
                        outdir)                    


                elif mode == "paftools_lift":

                    print "Running paftools liftover cross-species step " + species + " " + other_species                    

                    paftools_lift_cross_species_step(
                        cross_species_event_dicts,
                        other_species,
                        species,
                        paf_paths_dict,
                        species_dict,
                        paftools_path,
                        outdir)

                elif mode == "ioe":

                    print "Running ioe cross-species step " + species + " " + other_species                    

                    ioe_files_cross_species_step(
                            cross_species_event_dicts,
                            other_species,
                            species,
                            ioe_paths_dict,
                            event_dicts)

                elif mode == "exon_count":

                    print "Running exon count cross-species step " + species + " " + other_species                    

                    exon_count_cross_species_step(
                        cross_species_event_dicts,
                        other_species,
                        species,
                        ioe_paths_dict,
                        species_dict,
                        event_dicts)

                print "Completing event dict " + species + " " + other_species

                splice_lib.complete_event_dict(
                    cross_species_event_dicts[other_species][species],
                    inform_using_ri_events = False,
                    suppress_unique_edges = True)

                ##check for non-zero exon lists
                ##check for preservation of exon counts
                ##check for preservation of event type

                print "Checking cross mapping integrity " + species + " " + other_species

                check_mapped_event_integrity(
                    other_species, 
                    species, 
                    cross_species_event_dicts[other_species][species], 
                    event_dicts[species], 
                    mapping_failure_log,
                    id_translator)

                if mode == "minimap":

                    print "Checking for minimap multiply mapped events " + species + " " + other_species                                        

                    assess_multiply_mapped_events(
                            id_translator,
                            cross_species_event_dicts[other_species][species],
                            species,
                            other_species,
                            multimapper_log)                    

                #splice_lib.check_integrity(cross_species_event_dicts[other_species][species])

    mapping_failure_log.close()

    if mode == "minimap":
    
        multimapper_log.close()

    print "Merging event dicts by species"                  
            
    merged_event_dict = merge_event_dicts_by_species(
        event_dicts, 
        cross_species_event_dicts, 
        species_list, 
        outdir)

    print "Completing merged dict"                  

    complete_merged_dict(
        merged_event_dict, 
        species_dict, 
        add_gene_symbols, 
        transcript_gtfs)

    print "Renaming merged events"                  

    renamed_merged_event_dict = rename_merged_events(
        merged_event_dict, 
        species_list)

    print "Writing ambiguous mapping log"                  

    native_events = ambiguous_mapper_log_write(
        renamed_merged_event_dict,
        outdir)

    print "Writing native events log"                  

    write_native_events(
        native_events,
        outdir)

    print "Writing output"                  

    write_output(
        renamed_merged_event_dict,
        outdir,
        species_list)

    del cross_species_event_dicts
    del merged_event_dict
    del renamed_merged_event_dict



def parse_arguments(input_args):

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--event_file_prefixes", 
        type = str, 
        help = ("Comma-separated list of prefixes " +
                "to event gtfs, ioes, and fasta files. " + 
                "e.g. human_splice_lib_events,chimpanzee_splice_lib_events " + 
                "which will allow the program to find " + 
                "human_splice_lib_events.gtf, human_splice_lib_events.fa, " + 
                "chimpanzee_splice_lib_events.gtf and " + 
                "chimpanzee_splice_lib_events.fa from " + 
                "input directory"), 
        required = True)

    parser.add_argument(
        "--species_list", 
        type = str, 
        help = ("Comma-separated list of species names. " + 
               "Expected to be in same order as " + 
               "corresponding prefixes supplied in " + 
               "--event_file_prefixes argument"), 
        required = True)

    parser.add_argument(
        "--indir", 
        type = str, 
        help = ("Path to input directory in which " + 
                "event gtf and fasta files can be " + 
                "found (.fa extension expected for " + 
                "fasta files) with prefixes matching " + 
                "those supplied to the " + 
                "--event_file_prefixes argument"), 
        required = True)

    parser.add_argument(
        "--genome_fasta_paths", 
        type = str, 
        help = "Comma-separated paths to genome fasta " + 
               "files.  In same order as '--species_list' " + 
               "argument")

    parser.add_argument(
        "--chain_paths", 
        type = str, 
        help = ("Comma-separated paths to chain files. " + 
               "For each pair of species, a chain file " + 
               "should be supplied for both directions. " + 
               "So, for n-species, n(n-1) chain files " + 
               "should be provided.  These should be sorted " + 
               "first by the 'from' species, and next by the " + 
               "'to' species, with each sorting in same order " + 
               "as '--species_list' argument.  E.g. for three " + 
               "species such as human,chimp,orang, the chain " + 
               "paths should be ordered as follows: " + 
               "/path/to/human_to_chimp,/path/to/human_to_orang," + 
               "/path/to/chimp_to_human,/path/to/chimp_to_orang," + 
               "/path/to/orang_to_human,/path/to/orang_to_chimp"))

    parser.add_argument(
        "--paf_paths",
        type = str,
        help = ("Comma-separated paths to paf files. " + 
               "For each pair of species, a paf file " + 
               "should be supplied for both directions. " + 
               "So, for n-species, n(n-1) chain files " + 
               "should be provided.  These should be sorted " + 
               "first by the 'from' species, and next by the " + 
               "'to' species, with each sorting in same order " + 
               "as '--species_list' argument.  E.g. for three " + 
               "species such as human,chimp,orang, the paf " + 
               "paths should be ordered as follows: " + 
               "/path/to/human_to_chimp,/path/to/human_to_orang," + 
               "/path/to/chimp_to_human,/path/to/chimp_to_orang," + 
               "/path/to/orang_to_human,/path/to/orang_to_chimp"))

    parser.add_argument(
        "--minimap_outdir", 
        type = str, 
        help = ("Path to directory for minimap output. " + 
                "Providing this argument also specifies " + 
                "that you want to run minimap2."))

    parser.add_argument(
        "--liftover_outdir", 
        type = str, 
        help = ("Path to directory for liftover output. " + 
               "Providing this argument also specifies " + 
               "that you want to run liftOver"))

    parser.add_argument(
        "--paftools_lift_outdir",
        type = str,
        help = ("Path to directory for paftools " + 
                "liftover output. Providing this argument " + 
                "indicates paftools lift run."))

    parser.add_argument(
        "--exon_count_outdir",
        type = str,
        help = ("Path to directory for exon count-based " + 
                "approach. As with ioe, this mode assumes that " + 
                "transcript ids are the same across species."))

    parser.add_argument(
        "--ioe_outdir",
        type = str,
        help = ("Path to directory for ioe mode output. " + 
                "ioe mode requires ioe file for each species " + 
                " to be passed to --ioe_file_paths. " + 
                "This mode assumes transcript ids are (at least partly) " + 
                "all species."))

    parser.add_argument(
        "--ioe_file_paths",
        type = str,
        help = ("Paths to ioe files as comma-separated list " + 
                "order must be same as passed to --species_list"))

    parser.add_argument(
        "--paftools_path",
        type = str,
        help = "Path to paftools executable",
        default = "paftools.js")    

    parser.add_argument(
        "--samtools_path", 
        type = str, 
        help = ("Path to samtools executable - default = " + 
               "'samtools'"), 
               default = "samtools")

    parser.add_argument(
        "--bedtools_path", 
        type = str, 
        help = "Path to bedtools executable - default = 'bedtools'", 
        default = "bedtools")

    parser.add_argument(
        "--bedToGenePred_path", 
        type = str, 
        help = ("Path to bedToGenePred executable - " + 
                "default = 'bedToGenePred'"), 
        default = "bedToGenePred")

    parser.add_argument(
        "--minimap2_path", 
        type = str, 
        help = ("Path to minimap2 executable - " + 
                "default = 'minimap2'"), 
        default = "minimap2")

    parser.add_argument(
        "--num_threads", 
        type = str, 
        help = "Number of threads for minimap2 to use - default = 1", 
        default = "1")

    parser.add_argument(
        "--genePredToGtf_path", 
        type = str, 
        help = ("Path to genePredToGtf executable - " + 
                "default = 'genePredToGtf'"), 
        default = "genePredToGtf")

    parser.add_argument(
        "--liftover_path", 
        type = str, 
        help = "Path to liftOver executable - default = 'liftOver'", 
        default = "liftOver")

    parser.add_argument(
        "--add_gene_symbols", 
        action = "store_true", 
        help = ("If set, the genes to which the events " + 
                "belong will be included.  Requires " + 
                "--gene_gtfs and can make use of " + 
                "--transcript_gtfs"))

    parser.add_argument(
        "--gene_gtfs", 
        type = str, 
        help = ("Comma-separated list of paths to " + 
                "annotation gtfs for each species, " + 
                "in the same order as --species_list. " + 
                "Required (but only used for) " + 
                "--add_gene_symbols or --exon_count_outdir"))

    parser.add_argument(
        "--transcript_gtfs", 
        type = str, 
        help = ("Comma-separated list of paths " + 
                "to transcript gtfs for each species " + 
                "in the same order as --species_list. " + 
                "Can be made use of by --add_gene_symbols " + 
                "to increase the number of events for " + 
                "which symbols can be identified"))

    args = parser.parse_args()


    if args.add_gene_symbols and not args.gene_gtfs:

        sys.exit("Passage of --add_gene_symbols " + 
                 "requires --gene_gtfs")


    ##require that either minimap_outdir OR liftover_outdir be set

    if (args.minimap_outdir is None and 
        args.liftover_outdir is None and 
        args.paftools_lift_outdir is None and 
        args.ioe_outdir is None and 
        args.exon_count_outdir is None):

        sys.exit("None of --minimap_outdir, " + 
                 "--liftover_outdir, --paftools_lift_outdir, --ioe_outdir " + 
                 "or --exon_count_outdir are set. At least one must be set")

    if args.minimap_outdir:

        if args.genome_fasta_paths is None:
            sys.exit("Use of minimap requires " + 
                     "genome fasta files.  Pass " + 
                     "comma-separated paths to " + 
                     "--genome_fasta_paths")

    if args.liftover_outdir:

        if args.chain_paths is None:
            sys.exit("Use of liftover requires chain " + 
                     "files.  Pass comma-separated paths " + 
                     "to --chain_paths")

    if args.paftools_lift_outdir:

        if args.paf_paths is None:
            sys.exit("Use of paftools liftover requires paf " + 
                     "files.  Pass comma-separated paths " + 
                     "to --paf_paths")   

    if args.ioe_outdir:

        if args.ioe_file_paths is None:
            sys.exit("Use of ioe requires ioe " + 
                     "files.  Pass comma-separated paths " + 
                     "to --ioe_file_paths")     

    if args.exon_count_outdir:

        if (args.gene_gtfs is None or 
            args.ioe_file_paths is None):

            sys.exit("Use of exon count requires gene_gtf " + 
                     "files.  Pass comma-separated paths " + 
                     "to --gene_gtfs")  


    return args


def main(args):


    if args.gene_gtfs:

        annotated_transcript_gtf_list = args.gene_gtfs.split(",")

        if args.transcript_gtfs:

            merged_stringtie_gtf_path_list = args.transcript_gtfs.split(",")



    event_dicts = {}
    species_dict = {}

    species_list = args.species_list.split(",")
    prefixes = args.event_file_prefixes.split(",")


    ## Make output directory

    if args.minimap_outdir:

        gm.run_bash_cmd(("mkdir -p " + args.minimap_outdir).split())
        genome_fasta_paths = args.genome_fasta_paths.split(",")


    if args.liftover_outdir:

        gm.run_bash_cmd(("mkdir -p " + args.liftover_outdir).split())
        chain_paths = args.chain_paths.split(",")
        ###place chain file paths in dict indexed first by from species, then by to species

        chain_paths_dict = {}


        for species in species_list:

            species_list_minus = [i for i in species_list if i != species]

            for other_species in species_list_minus:

                chain_paths_dict.setdefault(species, {}
                    ).update({other_species: chain_paths.pop(0)})


    if args.paftools_lift_outdir:

        gm.run_bash_cmd(("mkdir -p " + args.paftools_lift_outdir).split())
        paf_paths = args.paf_paths.split(",")
        ###place chain file paths in dict indexed first by from species, then by to species

        paf_paths_dict = {}


        for species in species_list:

            species_list_minus = [i for i in species_list if i != species]

            for other_species in species_list_minus:

                paf_paths_dict.setdefault(species, {}
                    ).update({other_species: paf_paths.pop(0)})


    if args.ioe_outdir:

        gm.run_bash_cmd(("mkdir -p " + args.ioe_outdir).split())
        ioe_file_paths = args.ioe_file_paths.split(",")
        ###place chain file paths in dict indexed first by from species, then by to species


    if args.exon_count_outdir:

        gm.run_bash_cmd(("mkdir -p " + args.exon_count_outdir).split())
        ioe_file_paths = args.ioe_file_paths.split(",")
        ###place chain file paths in dict indexed first by from species, then by to species


    if args.ioe_outdir or args.exon_count_outdir:

        ioe_paths_dict = {}

        for species in species_list:

            species_list_minus = [i for i in species_list if i != species]

            for other_species in species_list_minus:

                ioe_paths_dict.setdefault(species, ioe_file_paths.pop(0))



    for i, species in enumerate(species_list):

        species_dict[species] = {}

        if args.minimap_outdir:

            species_outdir = (
                args.minimap_outdir + 
                "/" + 
                species + 
                "_out/")

            species_dict[species]["minimap_outdir"] = species_outdir
            gm.run_bash_cmd(("mkdir -p " + species_outdir).split())
            species_dict[species]["genome_fasta"] = genome_fasta_paths[i]

        if args.liftover_outdir:

            species_outdir = (
                args.liftover_outdir + 
                "/" + 
                species + 
                "_out/")

            species_dict[species]["liftover_outdir"] = species_outdir
            gm.run_bash_cmd(("mkdir -p " + species_outdir).split())


        if args.gene_gtfs:

            print "Importing " + species + " gene GTF"

            species_dict[species]["annotated_transcript_gtf"] = annotated_transcript_gtf_list[i]

            species_dict[species]["annotated_transcript_dict"] = splice_lib.generate_standard_transcript_dict(
                annotated_transcript_gtf_list[i])

            splice_lib.sort_transcript_dict_exons(
                        species_dict[species]["annotated_transcript_dict"])
            splice_lib.add_junctions_to_transcript_dict(
                        species_dict[species]["annotated_transcript_dict"])

            if args.transcript_gtfs:

                print "Importing " + species + " transcript GTF"

                species_dict[species]["merged_stringtie_gtf_path"] = merged_stringtie_gtf_path_list[i]


        species_dict[species]["prefix"] = prefixes[i]

        species_dict[species]["native_event_gtf_path"] = (
            args.indir + 
            "/" + 
            species_dict[species]["prefix"] + 
            ".gtf")

        try: 
            with open(species_dict[species]["native_event_gtf_path"]) as file:

                pass

        except IOError:

            sys.exit(species_dict[species]["native_event_gtf_path"] + 
                     " not found. Event GTF required for all species. " + 
                     "Exiting . . . ")



        if args.minimap_outdir:

            species_dict[species]["event_fasta_path"] = (
                args.indir + 
                "/" + 
                species_dict[species]["prefix"] + ".fa")

            try:
                with open(species_dict[species]["event_fasta_path"]) as file:
                    pass
            except IOError:
                sys.exit(species_dict[species]["native_event_gtf_path"] + 
                         " not found. Event fasta required for all " + 
                         "species for minimap2-based approach.  Exiting . . . ")

        print "Importing " + species + " event GTF"

        event_dicts[species] = splice_lib.generate_standard_event_dict(
            species_dict[species]["native_event_gtf_path"])

        print "Completing " + species + " event dict"

        splice_lib.complete_event_dict(event_dicts[species],
                                        inform_using_ri_events = False,
                                        suppress_unique_edges = True)


    if args.minimap_outdir:

        print "Running minimap cross-species processes"

        run_cross_species_processes(
            species_dict = species_dict,
            species_list = species_list,
            event_dicts = event_dicts,
            outdir = args.minimap_outdir,
            num_threads = args.num_threads,
            minimap2_path = args.minimap2_path,
            samtools_path = args.samtools_path,
            bedToGenePred_path = args.bedToGenePred_path,
            genePredToGtf_path = args.genePredToGtf_path,
            bedtools_path = args.bedtools_path,
            add_gene_symbols = args.add_gene_symbols,
            transcript_gtfs = args.transcript_gtfs,
            mode = "minimap",
            chain_paths_dict = None,
            liftover_path = None,
            paf_paths_dict = None,
            paftools_path = None,
            ioe_paths_dict = None)
   

    if args.liftover_outdir:

        print "Running UCSC liftover cross-species processes"

        run_cross_species_processes(
            species_dict = species_dict,
            species_list = species_list,
            event_dicts = event_dicts,
            outdir = args.liftover_outdir,
            num_threads = None,
            minimap2_path = None,
            samtools_path = None,
            bedToGenePred_path = None,
            genePredToGtf_path = None,
            bedtools_path = args.bedtools_path,
            add_gene_symbols = args.add_gene_symbols,
            transcript_gtfs = args.transcript_gtfs,
            mode = "liftover",
            chain_paths_dict = chain_paths_dict,
            liftover_path = args.liftover_path,
            paf_paths_dict = None,
            paftools_path = None,
            ioe_paths_dict = None)

    if args.paftools_lift_outdir:

        print "Running paftools liftover cross-species processes"

        run_cross_species_processes(
            species_dict = species_dict,
            species_list = species_list,
            event_dicts = event_dicts,
            outdir = args.paftools_lift_outdir,
            num_threads = None,
            minimap2_path = None,
            samtools_path = None,
            bedToGenePred_path = None,
            genePredToGtf_path = None,
            bedtools_path = None,
            add_gene_symbols = args.add_gene_symbols,
            transcript_gtfs = args.transcript_gtfs,
            mode = "paftools_lift",
            chain_paths_dict = None,
            liftover_path = None,
            paf_paths_dict = paf_paths_dict,
            paftools_path = args.paftools_path,
            ioe_paths_dict = None)

    if args.ioe_outdir:

        print "Running ioe cross-species processes"      

        run_cross_species_processes(
            species_dict = species_dict,
            species_list = species_list,
            event_dicts = event_dicts,
            outdir = args.ioe_outdir,
            num_threads = None,
            minimap2_path = None,
            samtools_path = None,
            bedToGenePred_path = None,
            genePredToGtf_path = None,
            bedtools_path = None,
            add_gene_symbols = args.add_gene_symbols,
            transcript_gtfs = args.transcript_gtfs,
            mode = "ioe",
            chain_paths_dict = None,
            liftover_path = None,
            paf_paths_dict = None,
            paftools_path = None,
            ioe_paths_dict = ioe_paths_dict)        

    if args.exon_count_outdir:

        print "Running exon count cross-species processes"       

        run_cross_species_processes(
            species_dict = species_dict,
            species_list = species_list,
            event_dicts = event_dicts,
            outdir = args.exon_count_outdir,
            num_threads = None,
            minimap2_path = None,
            samtools_path = None,
            bedToGenePred_path = None,
            genePredToGtf_path = None,
            bedtools_path = None,
            add_gene_symbols = args.add_gene_symbols,
            transcript_gtfs = args.transcript_gtfs,
            mode = "exon_count",
            chain_paths_dict = None,
            liftover_path = None,
            paf_paths_dict = None,
            paftools_path = None,
            ioe_paths_dict = ioe_paths_dict)   




if __name__ == '__main__':

    output_args = parse_arguments(sys.argv[1:])
    main(output_args)
