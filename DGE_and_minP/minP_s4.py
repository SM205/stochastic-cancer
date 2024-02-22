from multiprocessing import Pool
import pandas as pd
import numpy as np
import random
from scipy.stats import mannwhitneyu

random.seed(205)
s4h = pd.read_csv("../../Differential_Expression/True_nonzero/tmm_s4_periphery.csv",index_col="Unnamed: 0")
s4t = pd.read_csv("../../Differential_Expression/True_nonzero/tmm_s4_tumour.csv",index_col="Unnamed: 0")
s4 = s4h.join(s4t)
perm = 10000
healthy_perm = [None]*perm
tumour_perm = [None]*perm
for i in range(perm):
    healthy_perm[i] = random.sample(list(s4.columns),len(s4h.columns))
    tumour_perm[i] = list(set(s4.columns) - set(healthy_perm[i]))
genes = list(s4.index)
alphabetical_p = {}
for gene in genes:
    alphabetical_p[gene]=mannwhitneyu(s4h.loc[gene],s4t.loc[gene])[1]
ordered_p = dict(sorted(alphabetical_p.items(), key=lambda item: item[1]))
ascending_genes = list(ordered_p.keys())
ascending_p = list(ordered_p.values())
log2fc = list(np.log2((s4t.mean(axis=1)+1)/(s4h.mean(axis=1)+1))[ascending_genes])


def minP(i):
    p = []
    h_temp = s4[healthy_perm[i]]
    t_temp = s4[tumour_perm[i]]
    for gene in ascending_genes:
        p.append(mannwhitneyu(h_temp.loc[gene],t_temp.loc[gene])[1])
    
    q = [None]*len(p)
    p.reverse()
    for j in range(len(p)):
        if j > 0:
            q[j] = min(q[j-1],p[j])
        else:
            q[j] = p[j]           
    q.reverse()
    count = [0]*(len(s4.index))
    for m in range(len(p)):
        if q[m] <= ascending_p[m]:
            count[m] = count[m] + 1
    return count

if __name__ == "__main__":
    with Pool() as pool:
        list_count = pool.map(minP, range(perm))
    count_array = np.row_stack(list_count)
    count = np.sum(count_array,axis=0).tolist()
    p = [None]*(len(s4.index))
    for i in range(len(s4.index)):
        p[i] = count[i]/perm
        if i > 0:
            p[i] = max(p[i-1],p[i])
    minP_result = pd.DataFrame({"gene_name":ascending_genes,"log 2 FC":log2fc,"orginial p":ascending_p,"adjusted p":p})
    minP_result.set_index("gene_name",inplace=True)
    minP_result.to_csv("s4_DE_permutation_wilcoxon.csv")