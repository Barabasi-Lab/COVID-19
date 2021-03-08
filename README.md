# README
3/1/2021

Code used for the analysis of the paper:

**Network Medicine Framework for Identifying Drug Repurposing Opportunities for COVID-19.**

*Deisy Morselli Gysi, Ítalo Do Valle, Marinka Zitnik, Asher Ameli, Xiao Gan, Onur Varol, Susan Dina Ghiassian, JJ Patten, Robert Davey, Joseph Loscalzo, and Albert-László Barabási.*

## Disease Module Definition

-   Code: `00_Intro.R`
-   Description: Identify the Largest Connected Component (LCC, Disease Module), calculates its significance, plot the network and the distribution of the random (Figure 2).

### Input

-   Protein-Protein Interaction (PPI): `DatasetS2.csv`
-   SARS-CoV-2 Targets: `DatasetS1.csv`

## Disease Separation

-   Code: `01_Disease_Separation.ipynb`
-   Description: Calculation of the network separation $S_{AB}$ between diseases and the COVID-19 module (data presented in Figures S1 and S2).

### Input

-   Protein-Protein Interaction (PPI): `DatasetS2.csv`
-   SARS-CoV-2 Targets: `DatasetS1.csv`
-   Disease-Associated Genes: `Guney2016_GenesDisease.tsv`

## Drug repurpusing Pipelines

### AI

-   Code: `02_DR_AI.py`
-   Description: Calculates the Drug Repurpusing Candidates for pipelines A1-A4.

### Diffusion

-   Codes: `03_1_DR_Diff_matrix_generator.py` and `03_2_DR_Diff_run_pipelines.py`
-   Description: Calculates the Drug Repurpusing Candidates for pipelines D1-D5.

#### Input

-   Protein-protein interaction network: The graph needs to be a single connected component. Nodes are entrez IDs and in the format of an integer. File: `DatasetS2.csv`
-   Seeds: Disease-associated genes in the format of gene symbols. Each gene must be in one row similar to the sample file located at `DatasetS1.csv`
-   drug2targets: A dictionary of drugs (keys) and the targets (values). The targets are shown with Entrez IDs and are in the format of an integer. A sample file is located at `DatasetS3.csv`
-   app_drugs are similar to above. The app_drugs are the known approved drugs for the disease of study. They will be used as true positives to calculate the performance of the models. If there are no app_drugs, an empty dictionary {} should be fed into the code.

#### Running the code [FIX: It used the CT04, but we need to make it to S10-Sheet2]

1.  Before running the code, the user must generate the pre-calculated matrices, DSD, JS, and KL. To do so, `03_1_DR_Diff_matrix_generator.py` must be run.

-   `03_1_DR_Diff_matrix_generator.py` takes as inputs the network and will generate the similarity matrices and the corresponding node to index conversion dictionary and dumps the results under `output/diffusion`. They will be used later during the ranking process.

2.  Once the pre-calculated matrices are generated, `03_2_DR_Diff_run_pipelines.py` should be run. Results will be printed under `output/diffusion/out`.

#### Output

It is a `.csv` file with the ranking values of each model and the drug score. There is also a ranking value (integer) that corresponds to each ranking score. The file contains the following columns.

-   drug ID: the keys of the drug2targets dictionary
-   APP-Drugs: a binary value indicating whether or not the drug is approved. If app_drugs given by the user is empty, all binary values will be marked 'no'.
-   Degree: the number of targets
-   DSD-min-Rank: an integer showing drug rank by DSD-min; the lower, the better
-   KL-med-Rank: an integer showing drug rank by KL-med; the lower, the better\
-   KL-min-Rank: an integer showing drug rank by KL-min; the lower, the better\
-   JS-med-Rank: an integer showing drug rank by JS-med; the lower, the better\
-   JS-min-Rank: an integer showing drug rank by JS-min; the lower, the better\
-   DSD-min: the score given by DSD-min to each drug; the lower, the better\
-   KL-med: the score given by KL-med to each drug; the lower, the better\
-   KL-min: the score given by KL-min to each drug; the lower, the better\
-   JS-med: the score given by JS-med to each drug; the lower, the better\
-   JS-min: the score given by JS-min to each drug; the lower, the better

### Proximity

-   Code: `04_DR_Proximity.py`

-   Description: Calculates the Drug Repurpusing Candidates for pipelines P1-P3.

#### Input files and parameter

1.  disease gene list: a list of disease genes in entrez ID
2.  drug target file: a tab-delimited file with first entry being the drug, and the following entries being entrez IDs of targets
3.  interactome file: tab-delimited node pairs in their entrez ID. The 3rd column in the sample interactome is not being used
4.  -nsims specifies the number of randomizations in the proximity algorithm

#### Running the code

`> python 04_DR_Proximity.py -i_disease data/DatasetS1.csv -i_drug data/drug_targets_test.txt -i_interactome data/DatasetS2.csv -nsims 1000`

-   Sample Output display: (Note the beginning (first four lines) of this output are imported from previous code and are not relevant, please pay no attention to them).

<!-- -->

    >>>
    GenRev not found, steiner wont work
    Import error: Negex. Using keyword matching instead
    Import error: Funcassociate. Make sure that funcassociate is in toolbox!
    DIAMOnD not found and thus will not be available!
    Done parsing 332 disease genes

    > Done parsing drug targets: read 1 total  drugs

    > done loading network:
    > network contains 18505 nodes and 327924 links
    number of samples: 1
    drug=DBxxxxx, disease=DatasetS1.csv
    Finished in 0.019384 hours
    Current time: 
    18:51:56

#### Output files

The output of the code are .txt files in the `output/proximity` folder. If the \``output/proximity` folder is not present, it will automatically be created. Each `.txt` file is named with drug-disease input pair, e.g. `<Drug_ID>_<Disease_genes_input_filename>.txt` in the toy example. The output file consists of six tab-delimited entries. The entries are (respectively):

(1) drug,
(2) disease file,
(3) proximity distance, d,
(4) proximity Z-score, z,
(5) mean of proximity distance d in random simulation,
(6) standard deviation of proximity distance d in random simulation

#### Notes [Fix names, add P3]

1.  In the example `drug_targets_test.txt` file, we provide a made-up toy example drug for testing and demonstrating the code. The drug target data for each pipeline can be found in `drug_targets_P1.txt` `drug_targets_P2.txt`, and `Dataset_S4.txt`.

2.  computations with many drug-disease pairs are heavy. Parallel computing is recommended. This code has implementation of parallel computing, but it is not used in the demo above. Please check out the source code in `04_DR_Proximity.py`.

3.  The randomization of proximity algorithm is, by its mathematical nature, stochastic, i.e. re-running it does NOT generate exactly the same result for its Z-Score. Details can be found in ref 1. Guney(2016) paper.

The code here uses the code published in：

1.  Guney, E., Menche, J., Vidal, M. et al. Network-based in silico drug efficacy screening. Nat Commun 7, 10331 (2016). <https://doi.org/10.1038/ncomms10331>

2.  Menche J, Sharma A, Kitsak M, et al. Disease networks. Uncovering disease-disease relationships through the incomplete interactome. Science. 2015;347(6224):1257601. <doi:10.1126/science.1257601>


After the proximity code is run, you have to combine all `txt` files, they can be easily combined by `cat *.txt >> output.txt `.
