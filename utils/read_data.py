#######################################
### Project: Network Medicine Framework for Identifying Drug Repurposing Opportunities for COVID-19.
### Description: Read Data in for D1-D5
### Author: Asher Ameli
### email: asher dot ameli at gmail dot com 
### date: 1st March 2021
#######################################

import os

def connected_component_subgraphs(G, copy=True):
    import networkx as nx
        
    ## this function was removed from latest versions of networkx!!
    for c in nx.connected_components(G):
        if copy:
            yield G.subgraph(c).copy()
        else:
            yield G.subgraph(c)
            
            
def PPI():

    import networkx as nx

    path_cur = os.getcwd()
    #path_par = os.path.abspath(os.path.join(path_cur, os.pardir))
    path_par = path_cur
    file = open(path_par+'/data/DatasetS2.csv', 'r')
    rows = file.read().splitlines()[1::]
    file.close()
    G_o = nx.Graph()

    for row in rows:

        n = row.strip().split(',')[0:2]
        G_o.add_node(int(n[0]))
        G_o.add_node(int(n[1]))
        G_o.add_edge(int(n[0]), int(n[1]))

    clusters = sorted(list(connected_component_subgraphs(G_o)), key=len, reverse=True)
    G = clusters[0]
    G.remove_edges_from(nx.selfloop_edges(G))

    return G


def seed_gordon():

    import pandas as pd
    path_cur = os.getcwd()
    path_par = path_cur #os.path.abspath(os.path.join(path_cur, os.pardir))
    genes = list(pd.read_csv(path_par+'/data/DatasetS1.csv')['Symbol'].str.strip().values)

    return set(convert2entrez(genes))


def entrez_sym_map():

    path_cur = os.getcwd()
    #path_par = os.path.abspath(os.path.join(path_cur, os.pardir))
    path_par = path_cur
    file = open(path_par+'/data/interactome_2019_merged_protAnnots.csv', 'r')
    rows = file.read().splitlines()[1::]
    file.close()

    gene2entrez=dict()
    entrez2gene=dict()

    for row in rows:
        n = row.strip().split(',')
        entrez = int(n[1])
        sym = n[2]
        gene2entrez[sym] = entrez
        entrez2gene[entrez] = sym

    return gene2entrez, entrez2gene


def convert2entrez(symbols):

    gene2entrez, entrez2gene = entrez_sym_map()
    seed_list = []
    for sym in symbols:
        try:
            entrez = gene2entrez[sym]
        except KeyError:
            continue
        seed_list.append(entrez)

    return set(seed_list)


def drugbank():

    import pandas as pd
    path_cur = os.getcwd()
    #path_par = os.path.abspath(os.path.join(path_cur, os.pardir))
    path_par =  path_cur
    df = pd.read_csv(path_par + '/data/DatasetS3.csv')
    data= zip(df['ID'].values, df['entrez_id'].values)
    data_names = zip(df['ID'].values, df['Name'].values)
    drug2target={}
    for drug, entrez in data:
        try:
            drug2target.setdefault(drug, set()).add(int(entrez))
        except ValueError: continue

    drug2name = {}
    for drug, name in data_names:
        try:
            drug2name.setdefault(drug, set()).add(str(name).lower())
        except ValueError:
            continue

    results = {'drug2target': drug2target, 'drug2name':drug2name}

    return results


def clinical_trials():

    import pandas as pd
    path_cur = os.getcwd()
    #path_par = os.path.abspath(os.path.join(path_cur, os.pardir))
    path_par =  path_cur
    df = pd.read_excel(path_par+'/data/DatasetS10.xlsx',sheet_name = "DB")
    drugs = set(df['ID'].values)
    db_results = drugbank()
    drug2targets = db_results['drug2target']

    drug2targets_maped = {}
    for drug in drugs:
        drug2targets_maped[drug]= drug2targets[drug]

    return drug2targets_maped



