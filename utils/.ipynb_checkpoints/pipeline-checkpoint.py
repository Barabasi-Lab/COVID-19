class NetMeasure:

    def __init__(self, G, seed_ids, drugs_to_targets, app_drugs_to_targets):

        net_nodes = set(G.nodes())
        drugs_list = list(drugs_to_targets.keys()) + list(app_drugs_to_targets.keys())

        targets_all = [frozenset(net_nodes&drugs_to_targets[drug]) for drug in drugs_to_targets.keys()] +\
                      [frozenset(net_nodes&app_drugs_to_targets[drug]) for drug in app_drugs_to_targets.keys()]

        self.drugs_all = drugs_list
        self.targets_all = targets_all
        self.app_drugs = [drug for drug in app_drugs_to_targets.keys()]

        targets_set_unique = set()
        for target in targets_all: targets_set_unique.add(target)

        self.net = G
        self.seeds = seed_ids & set(G.nodes())
        self.degree = list(map(lambda target: len(target), self.targets_all))

        self.dsd_min = None
        self.dsd_med = None
        self.kl_min = None
        self.kl_med = None
        self.js_min = None
        self. js_med = None

    def dsd(self):

        from utils.func import dsd_rank
        drug_2_prob_median, drug_2_prob_min = dsd_rank(self.seeds, self.targets_all)

        self.dsd_med = list(map(lambda drug: drug_2_prob_median[drug], self.targets_all))
        self.dsd_min = list(map(lambda drug: drug_2_prob_min[drug], self.targets_all))

        return self

    def kl(self):
        from utils.func import kl_rank
        drug_2_prob_median, drug_2_prob_min = kl_rank(self.seeds, self.targets_all)

        self.kl_med = list(map(lambda drug: drug_2_prob_median[drug], self.targets_all))
        self.kl_min = list(map(lambda drug: drug_2_prob_min[drug], self.targets_all))

        return self

    def js(self):

        from utils.func import js_rank
        drug_2_prob_median, drug_2_prob_min = js_rank(self.seeds, self.targets_all)

        self.js_med = list(map(lambda drug: drug_2_prob_median[drug], self.targets_all))
        self.js_min = list(map(lambda drug: drug_2_prob_min[drug], self.targets_all))

        return self

    def make_df(self):

        import pandas as pd
        import numpy as np

        dsd_res = self.dsd()
        kl_res = self.kl()
        js_res = self.js()

        data = {'Drug': self.drugs_all,
                'Degree': self.degree,'Targets': self.targets_all,
                'DSD-min': dsd_res.dsd_min,
                'KL-med': kl_res.kl_med, 'KL-min': kl_res.kl_min,
                'JS-med': js_res.js_med, 'JS-min': js_res.js_min}

        df = pd.DataFrame(data)

        models = ['DSD-min', 'KL-med', 'KL-min', 'JS-med', 'JS-min']

        for m in models:
            df[m + '-Rank'] = df[m].rank(method='dense', ascending=True)
            df[m + '-Percentage'] = df[m].rank(method='dense', pct=True, ascending=True)

        df['APP-Drugs'] = np.where(df['Drug'].isin(self.app_drugs), 'yes', 'no')
        df = df.drop_duplicates(subset='Drug', keep="first")
        df_final = df[['Drug', 'APP-Drugs', 'Degree', 'Targets'] + [m + '-Rank' for m in models] + [m for m in models]]

        return df_final
