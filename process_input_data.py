import os
import numpy as np
import json
import pandas as pd

def build_model_config(eta, tau, pathway_dict, compounds, data_observed_weights,
                       a_init,
                       landa, gamma, mu, observed_weight_vector,
                       serialize=False, json_filename=None):
    """
    Builds a model configuration and save to json file if appropriate

    Parameters
    ----------
    eta: 2D array, mapping from pathway to compounds
    tau: 2D array, mapping from compounds to masses
    pathway_dictionary: dictionary of 'pathway' and 'compounds'
    ....

    Returns
    -------
    model_config : dictionary
         dictionary containing all data and initialization needed to run PUMA
    """
    model_config = {
        "landa"   : landa,
        "gamma"   : gamma,
        "mu": mu,
        "a_init" : a_init.tolist(),
        "eta": eta.tolist(),
        "tau": tau.tolist(),
        "pathway_dict": pathway_dict,
        "compounds": list(compounds),
        "data_observed_weights": data_observed_weights,
        "observed_weight_vector": observed_weight_vector.tolist()
    }

    if serialize:
        with open(json_filename, 'w') as outfile:
            json.dump(model_config, outfile)

    return model_config


def load_json():
    """
    Loads the json file
    Returns
    -------
    model_config:
        data structure holding all info required to run PUMA
    """

    json_filename = "modelConfig.json"

    if 'PUMA_OUTPUT_DATA' in os.environ:
        outdata_dir = os.environ['PUMA_OUTPUT_DATA']
        json_file = os.path.join(outdata_dir, json_filename)
        if not (os.path.isfile(json_file)):
            print ("unknown json file")
            exit()
    else:
        print("PUMA Error: cannot find output dir")

    with open(json_file) as data_file:
        model_config = json.load(data_file)

    model_config["a_init"] = np.asarray(model_config["a_init"])
    model_config["tau"] = np.asarray(model_config["tau"])
    model_config["eta"] = np.asarray(model_config["eta"])
    model_config["observed_weight_vector"] = np.asarray(model_config["observed_weight_vector"])
    return model_config


def initialize_pathway_activity(eta):
    """
    Initializes pathway activities; default and necessary: all pathways are initialized to true

    Parameters
    ----------
    eta: 2D array, mapping from pathway to compounds

    Returns
    -------
    a_init : 1D array
        initialization of activity for all pathways
    """
    a_init = np.ones(eta.shape[0])
    return a_init


def create_compound_to_masses_tau_matrix(met_dict, compounds):
    """
    Create an f by k matrix mapping of compounds to masses

    Parameters
    ----------
    met_dict: ...
    compounds:...

    Returns
    -------
    tau : 2D array
        Binary mapping of compounds to Masses.  Compounds are sorted alphabetically. Masses are sorted in ascending
        order
    mass_values: list
        Sorted list (ascending order) of all masses in the met_dict
    """
    assert(len(compounds) == len(met_dict))

    mass_values = set([x[0] for x in met_dict.values()])
    mass_values = sorted(mass_values)

    f = len(compounds)
    k = len(mass_values)
    tau = np.zeros((f, k))

    for f in compounds:
        # go over mass values
        # assign tau to 1, for each compound that has a mass k
        the_mass = met_dict[f][0]
        k = mass_values.index(the_mass)
        tau[compounds.index(f)][k] = 1
    return tau, mass_values


def create_observed_weight_vector(data_observed_weights, mass_values, ppm =15):
    """
    Create an 1 by k matrix indicating if a mass was observed.

    Parameters
    ----------
    data_observed_weights:
        input data of all observed masses
    mass_values:
    tau:

    Returns
    -------
    observed_weight_vector: 1D array, 1xk
         binary vector if a particular mass is observed in the data_observed_weights
    """

    # the observed weight vector should be of the same width as the tau matrix k entries.
    # The observed weight vector is organized from low mass to high, with 1 entries for observed weight

    k = len(mass_values)
    observed_weight_vector = np.zeros((1, k))
    keep_observed_masses = []
    n_observed_mets_in_model = 0
    for x in mass_values:
        for w in data_observed_weights:
        # an observed mass may map to one or more compound
        # if metabolite mass falls within w+15ppm or w-15ppm, then map w to met
            fifteenppm = w * ppm / 1000000.
            found =  x >= w - fifteenppm and x <= w + fifteenppm
            if found:
                keep_observed_masses.append(x)
                observed_weight_vector[0][mass_values.index(x)] = 1
                n_observed_mets_in_model+=1
    return observed_weight_vector[0]


def create_pathways_metabolite_eta_matrix(pathway_dict, compounds):
    """
    Create eta: a map from pathways to metabolites
    Parameters
    ----------
    pathway_dict:
        pathway_dict: dictionary of {'pathway': ['p1', 'p2', .., 'compounds': ['C00001 C00002 C00003', 'C000x C000x']}
    compounds:
         sorted unique list of compounds from the pathways  that have mass measurements

    Returns
    -------
    eta:
    """
    p = len(pathway_dict['pathway'])
    f = len(compounds)
    eta = np.zeros((p, f))
    for i in range(p):
        mets_in_pathway =pathway_dict['compounds'][i].split()
        for m in mets_in_pathway:
            eta[i][compounds.index(m)] = 1
    return eta



def process_input_files(pathway_file, met_file, obs_file, serialize=False, json_filename = None):
    """
        Process input files to create data structures to run PUMA
    :param pathway_file:
        # Format:  each line has pathway name, and metabolites in the pathway
        # first row is a header: pathway, compounds
        # cge00010,C00036 C15973 C00111 C00103 C00022 C00118
        #  every compound to be used in PUMA should belong to 1 or more pathways based on reading this file.
    :param met_file:
        # Format: each line has a compound ID followed by mass as specified by network model
        # C00015,404.0022
        # or 'C00015', 404.0022 # does not work yet
        # first row is a header: compound, mass
    :param obs_file:
        # merged_observations_file
        # if multiple sample files, this is the merged observations;  there is a likely compound identified by other tools,
        # and a confidence score.
        # row zero is a header
        # Column 0 will be the mass; last column will be the most likely match found using other tools (but not used
        #   when processing the input files as we do here. here we only care about the precursor mass   )
        # Format: Precursor Mass	Retention Time	Observation Confidence	Average Observation Score	Compound
        # 176.0565	12.68	0.75	0.75	C20844
     :param serialize:
        #if true create model_config
     :param json_filename:

    Returns
    -------
    eta: 2D numpry array
      Binary mapping of pathways to compounds. Compounds are sorted alphabetically.
    tau : 2D array
       Binary mapping of compounds to masses.  Compounds are sorted alphabetically. Masses are sorted in ascending
    pathway_dict: dict with 'pathway' and 'compounds' keys
        values are lists of pathways, and corresponding lists of compounds
    compounds: list
        sorted listing of compounds
    data_observed_weights:  list
         sorted listing of all observed masses
    observed_weight_vector: 1D array
        Binary vector indicating if a particular mass is observed
    """

    # process pathway_file
    pathway_df = pd.read_csv (pathway_file)
    pathway_dict = pathway_df.to_dict('list')
    compounds_listing = list(pathway_dict['compounds'])

    #find unique compound listing
    compounds_list = sorted(set([x for item in compounds_listing for x in item.split()]))

    # process met_file
    met_df = pd.read_csv (met_file)
    met_dict = met_df.set_index('compound').T.to_dict('list')
    mets_with_masses = set([x for x in met_dict.keys()])

    # Check consistency of pathway and met files:
    compounds_listing_in_pathways = set(compounds_list)
    # Are there compounds without masses?  In compounds_listing_in_pathways but not in mets_with_masses
    compounds_with_no_masses = compounds_listing_in_pathways - mets_with_masses

    if (len(compounds_with_no_masses)>0 ):

        # print ("PUMA Warning 1: Compound in pathway does not have a matching mass. Please specify a mass for each "
        #        "compound"
        #        "For now %d compounds will be removed:" % len (compounds_with_no_masses))
        # print('compounds_with_no_masses:', list(compounds_with_no_masses))

        pathway_dict, compounds_listing = \
            return_sanitized_pathways_and_compounds_listing(pathway_dict, compounds_with_no_masses)

        compounds_list = sorted(set([x for item in compounds_listing for x in item.split()]))

        # are there still any issues that cannot be resolved?
        if(len(set(compounds_list) - mets_with_masses)>0):
            print("PUMA ERROR 6: still have a problem removing compounds that have no masses")
            exit()

    # Are there mass compounds listed in the masses file without a mention in the pathways file?
    compounds_not_in_pathways = mets_with_masses-set(compounds_list)

    if (len(compounds_not_in_pathways) > 0):
        # print("PUMA Warning 2: There exists  mass compounds in %s,\n but compounds not in pathways."
        #        "Please remove them from mass listing,"
        #        "or they will be ignored" % met_file)
        # print("%d compounds are not in pathways; however the masses are given" %len(compounds_not_in_pathways))
        #delete them from met dictionary
        for k in compounds_not_in_pathways:
            del met_dict[k]

    # Create mapping from pathways to metabolites
    eta = create_pathways_metabolite_eta_matrix(pathway_dict, compounds_list)
    # create mapping from metabolites to masses
    tau, unique_masses = create_compound_to_masses_tau_matrix(met_dict, compounds_list)

    print("pathways_in_model:", pathway_dict['pathway'])
    print("number_pathways_in_model:", len(pathway_dict['pathway']))
    print("metabolites_in_model:", compounds_list)
    print("number_metabolites_in_model:", len(compounds_list))

    # run santify check over matrices
    sanity_check(eta, tau)

    # create observed_weight_vector
    obs_df = pd.read_csv (obs_file)
    data_observed_weights= obs_df['mz'].values.tolist()
    data_observed_weights= sorted(data_observed_weights)

    # Note: we no longer need tau to pass to this function
    observed_weight_vector = create_observed_weight_vector(data_observed_weights, unique_masses)
    n_observed_masses = len([observed_mass for observed_mass in observed_weight_vector if observed_mass])

    return eta, tau, pathway_dict, compounds_list, data_observed_weights, observed_weight_vector


def return_sanitized_pathways_and_compounds_listing(pathway_dict, compounds_with_no_masses):
    """
    remvoing pathways which don't have any metabolites with corresponding masses
    Parameters
    ----------
    pathway_dict:
        pathway_dict: dictionary of {'pathway': ['p1', 'p2', .., 'compounds': ['C00001 C00002 C00003', 'C000x C000x']}
    compounds:
         sorted unique list of compounds from the pathways  that have mass measurements
    Returns
    -------
    santized pathways, compounds_listing

    """
    compounds_listing = []
    for item in pathway_dict['compounds']:
        new_value = set([x for x in item.split()])
        item = new_value - compounds_with_no_masses
        compounds_listing.append(' '.join(str(e) for e in item))

    pathway_dict['compounds'] = compounds_listing

    sanitized_compounds_listing = []
    to_remove_indices = []
    for idx, each in enumerate(compounds_listing):
        if each == "":
            # print("pathway %d has no metabolites with corresponding mass, so it will be removed" % idx, '\n')
            to_remove_indices.append(idx)
        else:
            sanitized_compounds_listing.append(each)

    compounds_listing = sanitized_compounds_listing
    pathway_dict['compounds'] = compounds_listing

    santized_pathway_list = [pathway for index, pathway in enumerate(pathway_dict["pathway"])
                             if index not in to_remove_indices]
    pathway_dict["pathway"] = santized_pathway_list
    return pathway_dict, compounds_listing


def sanity_check(eta, tau):

    # sanity check of data
    pm_row_sum = np.sum(eta, axis=0)
    pm_col_sum = np.sum(eta, axis=1)

    if np.any(np.min(pm_col_sum) <= 0):
        raise Exception('Column sum of the pathway-metabolite matrix has zero entries. ' +
                        'It means that some pathway doesn\'t generate metabolites.')
    if np.any(np.min(pm_row_sum) <= 0):
        raise Exception('Row sum of the pathway-metabolite matrix has zero entries. ' +
                        'It means that some metabolite is not generated by any pathways.')

    # sanity check of data
    mw_row_sum = np.sum(tau, axis=0)
    mw_col_sum = np.sum(tau, axis=1)

    if np.any(mw_col_sum != 1):
        raise Exception('Column sum of the metabolite-mass matrix has entries not equal to 1. ' +
                        'It means that some metabolite has no mass or more than one mass values.')

    if np.any(np.min(mw_row_sum) <= 0):
        raise Exception('Row sum of the metabolite-mass matrix has zero entries. ' +
                        'It means that some mass is not from any metabolite.')



def build_config(serialize=False, landa=0.5, gamma=0.9, mu=0.5,
                 pathway_filename = "pathways.csv", met_filename = "masses.csv",
                 observations_filename = "observed_data.csv", json_filename = "modelConfig.json"):

    """
    Builds configuration for running PUMA using input files
    :param serialize:
        if true, write out a json file with model configuration
    :param landa:
        probability of a pathway being active.  Currently same value for all pathways
    :param gamma:
              
    :param mu:
    :param pathway_filename:
    :param met_filename:
    :param observations_filename:
    :param json_filename:
    Rreturns
    --------
    model_config:
        data structure holding all info required to run PUMA
    """

    # setup path for each input file
    if 'PUMA_INPUT_DATA' in os.environ:
        indata_dir= os.environ['PUMA_INPUT_DATA']
        pathway_file = os.path.join(indata_dir, pathway_filename)
        met_file = os.path.join(indata_dir, met_filename)
        obs_file = os.path.join(indata_dir, observations_filename)
        # check that all input files exist:
        if not (os.path.isfile(pathway_file)):
            print ("PUMA Error 1: unknown pathways file")
        if not (os.path.isfile(met_file)):
            print ("PUMA Error 2: unknown mets file")
        if not (os.path.isfile(obs_file)):
            print ("PUMA Error 3: unknown observation file")
        if not ((os.path.isfile(pathway_file)) and (os.path.isfile(met_file)) and os.path.isfile(obs_file)):
            exit()
    else:
        print ("PUMA Error 4: cannot find input data dir in environment")

    #  setup path for output files
    if serialize and 'PUMA_OUTPUT_DATA' in os.environ: # path exists
        outdata_dir= os.environ['PUMA_OUTPUT_DATA']
        json_file = os.path.join(outdata_dir, json_filename)
        if not os.path.exists(outdata_dir):  # check directory to write json file exists
            print ("PUMA Error 5: unknown output directory")

    # process the input files
    eta, tau, pathway_dict, compounds, data_observed_weights, observed_weight_vector = process_input_files(pathway_file,
                                                                                                           met_file,
                                                                                                           obs_file)
    # initialize pathway activity
    a_init = initialize_pathway_activity(eta)

    # build a model configuration
    model_config = build_model_config(eta, tau, pathway_dict, compounds, data_observed_weights,
                                      a_init,
                                      landa, gamma, mu, observed_weight_vector,
                                      serialize, json_filename=json_file)
    return model_config


def main():

    model_config = load_json()


if __name__ == "__main__":
    main()
