#######################################
### Project: Network Medicine Framework for Identifying Drug Repurposing Opportunities for COVID-19.
### Description: Pipeline for Proximity: P1-P3
### Author: Xiao Gan
### email: jack dot xiao dot gan at gmail dot com 
### date: 1st March 2021
#######################################

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os
import re

import sys
sys.path.append('./utils/')

from guney_code import wrappers
from guney_code import network_utilities
from guney_code import network_utils
#import testdisease

import separation as tools
from multiprocessing import Pool


"""
Neccesary files:
guney_code(folder)
genes_annotated_0311.csv
separation.py
disease genes file: interactome.tsv
disease genes file: COVID19_Human_Targets.csv
drug target file: drug_targets_test.txt

"""


def convert(list1,entry_from='GeneID',entry_to='Symbol' ,dataset = 'data/interactome_2019_merged_protAnnots.csv'):
    # convert a list of entry1 to another list of entry2
    df = pd.read_csv(dataset)
    listx1 = list(df[entry_from].astype(str))
    listx2 = list(df[entry_to].astype(str))
    result_list =[]
    for i in list1:
        try:
            result_list.append(listx2[listx1.index(i)])
        except:
            raise Exception('Error occured')
    return result_list

def parse_drug_target(file1, output_file=None):
    # read a file with 1st column of drug name and the other columns being target genes
    # returns a dictionary with keys as drug names and values as lists of its target genes
    # output_file into
    drug_targets={}
    if output_file != None:
        f = open(output_file,"a+")
    for line in open(file1):
        words = line.split("\t")
        set1=set()
        for word1 in words[1:]:
            if word1 != ('\n') and word1 != (''):
                set1.add(word1.strip('\n'))
        drug_targets[words[0]]=set1
        if output_file != None:
            f.write(str(words[0])+'\t'+ str(set1)+'\n')
    print (('\n> Done parsing drug targets: read %s total  drugs')%len(drug_targets))
    if output_file != None:
        f.close()
    return drug_targets

def single_proximity(sample):
    # --------------------------------------------------------
    #
    # LOADING NETWORK and DISEASE GENES
    #
    # --------------------------------------------------------
    drug_key = sample[0]
    disease_key = sample[1]
    drug_targets =sample[2]
    disease_genes =sample[3]
    network = sample [4]
    nsims = sample [5]

    disease_save = disease_key.split("/")[1]
    disease_save = re.sub(".csv$", "", disease_save)
        
    # multiple proximity loop. Call each key combination from the two dictionaries

    nodes_from =  set(drug_targets[drug_key]) & set(network.nodes())
    nodes_to =  set(disease_genes[disease_key]) & set(network.nodes())
    print (('drug=%s, disease=%s')%(drug_key, disease_key))

    if len(nodes_from) == 0 | len(nodes_to) == 0: # if no applicable target, stop
        return

    # computing proximity. Please set the parameters to proper values.
    d, z, (mean, sd) = wrappers.calculate_proximity(network, nodes_from, nodes_to, n_random=nsims, min_bin_size = 100, seed=None)
    #print (('d=%s, z=%s, (mean, sd)=%s')%(d, z, (mean, sd)))

    # write in file
    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, r'output/proximity/')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)


    filename= drug_key +'_'+ disease_save+'.txt'
    with open(os.path.join(final_directory, filename), 'w') as f:
        #f.write(drug_key +'\t'+ disease_key +'\td='+ str(d) +'\tz=' + str(z) +'\t(mean, sd)='+ str((mean, sd))+'\n')
        f.write(drug_key +'\t'+ disease_save +'\t'+ str(d) +'\t' + str(z) +'\t'+ str(mean) +'\t'+ str(sd)+'\n')

    return

import argparse
# Initialize parser
parser = argparse.ArgumentParser(description='Compute drug-disease proximity.')

# Adding optional argument
parser.add_argument("-i_drug", "--drug_targets_file", help = "input drug_target file")
parser.add_argument("-i_disease", "--disease_gene_file", help = "input disease gene file")
parser.add_argument("-i_interactome", "--network_file1", help = "input interactome file")
parser.add_argument("-nsims", "--nsims", help = "number of randomization")

# Read arguments from command line
args = parser.parse_args()
disease_gene_file = args.disease_gene_file
drug_targets_file = args.drug_targets_file
network_file1 = args.network_file1
nsims = int(args.nsims)
##print (args)

# get disease genes
df = pd.read_csv(disease_gene_file)
disease_gene_set = set(df['EntrezID'].astype(str))

print (('Done parsing %s disease genes')%len(disease_gene_set))

# get drug targets
##drug_targets_file = 'drug_targets_PPIlcc_0402.txt'
##drug_targets_file = 'drug_targets_test.txt'
drug_targets = parse_drug_target(drug_targets_file, output_file=None)

##network_file1 = 'interactome2019.txt'
G  = tools.read_network(network_file1)

# filter for nodes in lcc of G
components = nx.connected_components(G)
lcclist = sorted(list(components), key = len, reverse=True)
G1 = nx.subgraph(G, lcclist[0])

disease_gene_set = disease_gene_set & set(G1.nodes())
disease_genes=dict()
disease_genes[disease_gene_file] = disease_gene_set

samples=[]
for i in drug_targets.keys():
    if len(drug_targets[i]) >= 2:
        for j in disease_genes.keys():
            samples.append([i,j,drug_targets,disease_genes,G1,nsims])
print ('number of samples:',len(samples))

# set up parallel computing
ncpus=1

import time
import datetime
start_time = time.time()

s = time.time()
for sample in samples:
    single_proximity(sample)
##p = Pool(ncpus)
##res = p.map(single_proximity, samples)
##p.close()

e = time.time() - s

print ('Finished in %f hours'%(e/3600))

now = datetime.datetime.now()

print ("Current time: ")
print (now.strftime("%H:%M:%S"))

