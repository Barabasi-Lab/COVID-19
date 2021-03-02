#######################################
### Project: Network Medicine Framework for Identifying Drug Repurposing Opportunities for COVID-19.
### Description: Drug Repurpusing for D1-D5
### Author: Asher Ameli
### email: asher dot ameli at gmail dot com 
### date: 1st March 2021
#######################################

from utils.pipeline import NetMeasure
from utils.read_data import PPI, drugbank, seed_gordon, clinical_trials

# - - - - - - - - - - - - - - Network, Seeds, APP Targets - - - - - - - - - - -
drugbank_results = drugbank()
drug2targets = drugbank_results['drug2target']
seed_ids = seed_gordon()
app_drugs = clinical_trials()
G = PPI()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
network_status = NetMeasure(G, seed_ids, drug2targets, app_drugs)
df = network_status.make_df()
df.to_csv('output/diffusion/out/master-list.csv', index=False)


