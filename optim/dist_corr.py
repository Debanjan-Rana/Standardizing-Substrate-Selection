"""script for optimizing umap parameters based on the distance correlation metric D.
The parameter values can be specified as list. The script returns results as csv file for
all the different parameter settings and plots the results."""

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from correlation_metric import optimize_umap, save_boxplot, true_jaccard_dist
from pathlib import Path
import pandas as pd
from scipy import stats


ndims = [2]
neighbours = [15, 30]
min_dist = [0.1, 0.2]
cluster_size = 15
metric = 'jaccard'
folder = Path.cwd().joinpath("results", "test_rep")


smiles_df = pd.read_csv('/home/employee/rdebanja/pycharm_projects/SUBCLUST/Standardized subs '
                        'selection/data/final_DB/final_drugbank.csv')['SMILES']
mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_df]
fp_list = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 512) for mol in mol_list]
FP_array = np.array(fp_list)

# computing the jaccard distance between all fp pairs (original jaccard distance)
original_jaccard_dist = true_jaccard_dist(fp_list)
jacc_dist_vec = original_jaccard_dist[np.triu_indices(len(mol_list))]
norm_ori_jacc_dist = jacc_dist_vec / np.max(jacc_dist_vec)


data_store = np.zeros((len(ndims)*len(neighbours)*len(min_dist), 5))
count = 0
for ndim in ndims:
    for nb in neighbours:
        for mdist in min_dist:
            reduced_distance = optimize_umap(FP_array, ndim, nb, mdist, metric, cluster_size, folder)
            p_corr = stats.pearsonr(norm_ori_jacc_dist, reduced_distance)
            print(p_corr)
            title = f"Dim={ndim}_nb={nb}_m-dist={mdist}_score={p_corr[0]: .3f}"
            save_boxplot(norm_ori_jacc_dist, reduced_distance, title, folder, intervals=20)
            data_store[count, :] = [ndim, nb, mdist, p_corr[0], p_corr[1]]
            data_df = pd.DataFrame(data_store, columns=['dim', 'nbs', 'mdist', 'p_corr0', 'p_corr1'])
            data_df.to_csv('correlation_metric_analysis.csv', index=False)
            count += 1

print('asd')