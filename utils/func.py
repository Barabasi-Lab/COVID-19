import pickle
import numpy as np
import os
import warnings
warnings.filterwarnings("ignore")


def dsd_rank(seeds, set_drugs):

    path_cur = os.getcwd()
    #path_par = os.path.abspath(os.path.join(path_cur, os.pardir))
    path_par = path_cur
    
    node_2_id = pickle.load(open(path_par + '/output/diffusion/ppi_node_2_id.p', 'rb'))
    dsd_val = pickle.load(open(path_par + '/output/diffusion/PPI_DSD_100.p', 'rb'))

    node_2_prob_median = {}
    node_2_prob_min = {}
    # -------
    for drug in set_drugs:

        prob_all_med = []
        prob_all_min = []

        for target in drug:
            prob_drug = []
            for s in seeds:
                if target == s: continue
                prob = dsd_val[node_2_id[target], node_2_id[s]]
                prob_drug.append(prob)

            med_prob = np.median(prob_drug)
            min_prob = np.min(prob_drug)
            prob_all_med.append(med_prob)
            prob_all_min.append(min_prob)

        node_2_prob_median[drug] = np.mean(prob_all_med)
        node_2_prob_min[drug] = np.mean(prob_all_min)

    return node_2_prob_median, node_2_prob_min


def kl_rank(seeds, set_drugs):

    path_cur = os.getcwd()
    #path_par = os.path.abspath(os.path.join(path_cur, os.pardir))
    path_par = path_cur
    
    node_2_id = pickle.load(open(path_par + '/output/diffusion/ppi_node_2_id.p', 'rb'))
    kl_val = pickle.load(open(path_par + '/output/diffusion/PPI_KL_100.p', 'rb'))

    node_2_prob_median = {}
    node_2_prob_min = {}
    # -------
    for drug in set_drugs:

        prob_all_med = []
        prob_all_min = []

        for target in drug:
            prob_drug = []
            for s in seeds:
                if target == s: continue
                prob = kl_val[node_2_id[target], node_2_id[s]]
                prob_drug.append(prob)

            med_prob = np.median(prob_drug)
            min_prob = np.min(prob_drug)
            prob_all_med.append(med_prob)
            prob_all_min.append(min_prob)

        node_2_prob_median[drug] = np.mean(prob_all_med)
        node_2_prob_min[drug] = np.mean(prob_all_min)

    return node_2_prob_median, node_2_prob_min


def js_rank(seeds, set_drugs):

    path_cur = os.getcwd()
    #path_par = os.path.abspath(os.path.join(path_cur, os.pardir))
    path_par = path_cur

    node_2_id = pickle.load(open(path_par + '/output/diffusion/ppi_node_2_id.p', 'rb'))
    js_val = pickle.load(open(path_par + '/output/diffusion/PPI_JS_100.p', 'rb'))

    node_2_prob_median = {}
    node_2_prob_min = {}
    # -------
    for drug in set_drugs:

        prob_all_med = []
        prob_all_min = []

        for target in drug:
            prob_drug = []
            for s in seeds:
                if target == s: continue
                prob = js_val[node_2_id[target], node_2_id[s]]
                prob_drug.append(prob)

            med_prob = np.median(prob_drug)
            min_prob = np.min(prob_drug)
            prob_all_med.append(med_prob)
            prob_all_min.append(min_prob)

        node_2_prob_median[drug] = np.mean(prob_all_med)
        node_2_prob_min[drug] = np.mean(prob_all_min)

    return node_2_prob_median, node_2_prob_min
