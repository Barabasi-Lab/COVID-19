#######################################
### Project: Network Medicine Framework for Identifying Drug Repurposing Opportunities for COVID-19.
### Description: Matrix Generator for D1-D5
### Author: Asher Ameli
### email: asher dot ameli at gmail dot com 
### date: 1st March 2021
#######################################

from __future__ import division
import numpy as np
from scipy.spatial.distance import pdist, squareform
import networkx as nx
from scipy.stats import entropy
import pickle


def kl_divergence(p, q):
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))


def js_cal(p, q):
    m = (p + q) / 2
    return (entropy(p, m) + entropy(q, m)) / 2


def kl_divergence_matrix(rw_matrix_in):

    rw = np.array(rw_matrix_in)
    rw_sum_rows = rw.sum(axis=1)
    rw = rw / rw_sum_rows[:, np.newaxis]

    n = rw_matrix_in.shape[0]
    rw_kl = np.zeros((n, n))
    rw_js = np.zeros((n, n))

    for i in range(0, n):
        for j in range(0, n):
            rw_kl[i, j] = kl_divergence(rw[i, :], rw[j, :])
            rw_js[i, j] = js_cal(rw[i,:], rw[j,:])

    results = {'RWKL': rw_kl, 'RWJS': rw_js}
    return results


def rw_dsd_generator(adjacency_in, nrw):

    adjacency = np.asmatrix(adjacency_in)
    n = adjacency.shape[0]
    degree = adjacency.sum(axis=1)
    p = adjacency / degree
    c = np.eye(n)
    for i in range(nrw):
        c = np.dot(c, p) + np.eye(n)

    dsd = squareform(pdist(c, metric='cityblock'))
    results = {'RW': c, 'DSD': dsd}
    return results


if __name__ == "__main__":

    from utils.read_data import PPI
    G = PPI()
    # --------------------------------------------
    G_nodes = G.nodes()
    node_2_id = {}
    id2node = {}
    for ind, node in enumerate(G_nodes):
        node_2_id[node] = ind
        id2node[ind] = node
    pickle.dump(id2node, open('output/diffusion/ppi_id_2_node.p', 'wb'))
    pickle.dump(node_2_id, open('output/diffusion/ppi_node_2_id.p', 'wb'))
    # --------------------------------------------
    adjacency_ppi = np.array(nx.adjacency_matrix(G).todense())
    results_rw_dsd = rw_dsd_generator(adjacency_ppi, 100)
    rw_matrix = np.array(results_rw_dsd['RW'])
    dsd_matrix = np.array(results_rw_dsd['DSD'])

    pickle.dump(dsd_matrix, open('output/diffusion/PPI_DSD_100.p', 'wb'))
    # ------------------- KL & JS ------------
    results_kl_js = kl_divergence_matrix(rw_matrix)
    kl_matrix = np.array(results_kl_js['RWKL'])
    js_matrix = np.array(results_kl_js['RWJS'])

    pickle.dump(kl_matrix, open('output/diffusion/PPI_KL_100.p', 'wb'))
    pickle.dump(js_matrix, open('output/diffusion/PPI_JS_100.p', 'wb'))