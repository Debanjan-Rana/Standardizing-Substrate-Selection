"""script for computing the average silhuette score for different umap parameters and no. of clusters.
The parameter values can be specified as a list."""

from rdkit import Chem
from rdkit.Chem import AllChem
from umap import UMAP
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score

# variation of parameters to be analyzed
ndims = [2]
neighbours = [15, 30]
min_dist = [0.1, 0.2]
cluster_size = [15, 20]
metric = 'jaccard'
folder = Path.cwd().joinpath("results", "ECFP4_nb_mdist")

smiles_df = pd.read_csv('/home/employee/rdebanja/pycharm_projects/SUBCLUST/Standardized subs '
                        'selection/data/final_DB/final_drugbank.csv')['SMILES']

mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_df]
fp_list = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 512) for mol in mol_list]
FP_array = np.array(fp_list)


def save_umap_fig(UMAPspace, title, cls, folder):
    fig, ax = plt.subplots()
    plt.scatter(UMAPspace[:, 0], UMAPspace[:, 1], s=8, c=cls, alpha=0.5)
    ax.set_xlabel(r'UMAP1', fontsize=15)
    ax.set_ylabel(r'UMAP2', fontsize=15)
    ax.set_title(title)
    fig.tight_layout()
    filename = "umap_" + title + ".png"
    plt_location = folder.joinpath(filename)
    plt.savefig(plt_location, dpi=120, transparent=True, bbox_inches="tight")
    # plt.show()
    plt.clf()
    plt.close('all')


data_store = np.zeros((len(ndims) * len(neighbours) * len(min_dist) * len(cluster_size), 5))
count = 0

for ndim in ndims:
    for nb in neighbours:
        for mdist in min_dist:

            UMAPspace = UMAP(n_components=ndim, n_neighbors=nb, min_dist=mdist, metric=metric,
                             random_state=0).fit_transform(FP_array)

            UMAPspace = np.nan_to_num(UMAPspace)
            for cls_size in cluster_size:
                z = linkage(UMAPspace, method="ward")
                cls = fcluster(z, cls_size, criterion="maxclust")
                silhuette_score = silhouette_score(UMAPspace, cls)
                print(nb, mdist, cls_size, silhuette_score)

                title = f"Dim={ndim}_nb={nb}_m-dist={mdist}_cls_size={cls_size}_silhscore={silhuette_score}"
                save_umap_fig(UMAPspace, title, cls, folder)
                data_store[count, :] = [ndim, nb, mdist, cls_size, silhuette_score]
                data_df = pd.DataFrame(data_store, columns=['dim', 'nbs', 'mdist', 'cls_size', 'silhuette_score'])
                data_store_loc = folder.joinpath('silh_score_analysis.csv')
                data_df.to_csv(data_store_loc, index=False)
                count += 1
print('asd')
