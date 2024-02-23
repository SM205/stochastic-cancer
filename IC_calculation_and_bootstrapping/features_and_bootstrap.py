import pandas as pd
import numpy as np

def get_observed_variance(genes):
    return genes.sum(axis=0).var()

def get_pb_variance(genes):
    fx = genes.sum(axis=1) / genes.shape[1]
    return np.sum([(1-p)*p for p in fx])

def get_IC(obs_var, pb_var):
    return obs_var / pb_var

def bootstrap_IC(genes, bootstrap_n = 10000):
    bootstrap_populations = [genes.sample(frac=1.0, replace=True, axis=1) for _ in range(bootstrap_n)]
    pb_vars = np.array([get_pb_variance(x) for x in bootstrap_populations])
    obs_vars = np.array([get_observed_variance(x) for x in bootstrap_populations])
    ics = obs_vars / pb_vars

    return pd.Series({
        "mean": np.nanmean(ics),
        "lower": np.nanquantile(ics, q=0.025),
        "median": np.nanquantile(ics, q=0.5),
        "upper": np.nanquantile(ics, q=0.975),
    })

family_genes = pd.read_csv("familyset.csv",index_col='gene_id')
family_df = pd.read_csv("clean_pantherhuman.csv")
dich = pd.read_csv("dichotomised_genes.csv",index_col='gene_id')
l = []
fam_list = set(family_genes["family_id"])
for family in fam_list:
    dic = {}    
    genes = dich.reindex(index=family_genes[family_genes.family_id == family].index)
    exp_per_cell = genes.sum(axis=0)
    exp_per_gene = genes.sum(axis=1)
    dic['family_id'] = family
    dic['mean_gene_per_cell'] = exp_per_cell.mean()
    dic['std_gene_per_cell'] = exp_per_cell.std()
    dic['original_IC'] = get_IC(get_observed_variance(genes),get_pb_variance(genes))
    dic['no_genes_total'] = family_df.loc[family_df.family_id == family].shape[0]
    dic['no_genes_measured'] = len(dich.loc[genes.index].dropna())
    dic['non_zero_genes'] = (dich.loc[genes.index].dropna().sum(axis=1)>0).sum()
    dic["mean"],dic["lower"],dic["median"],dic["upper"] = bootstrap_IC(genes)
    if exp_per_gene.sum() > 0:
        mean_per_gene_expression = exp_per_gene/(exp_per_gene.sum())
        dic['shannon_index'] = -((mean_per_gene_expression*np.log(mean_per_gene_expression)).sum())
        dic['species_evenness'] = dic['shannon_index']/np.log(dic['non_zero_genes'])
    else:
        dic['shannon_index'] = 0
        dic['species_evenness'] = 0
    dic['family_name'] = family_df[family_df.family_id==family].iloc[0]["family_name"]
    l.append(dic)
results = pd.DataFrame(l)
results.set_index("family_id",inplace=True)
results.to_excel("Final_IC_bootstrap.xlsx")