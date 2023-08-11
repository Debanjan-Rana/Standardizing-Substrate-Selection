"""functions to be used for computing the distance coprrelation metric for optimizing umap parameters"""

from rdkit import DataStructs
from umap import UMAP
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score
import seaborn as sns
from scipy.spatial import distance


def optimize_umap(FP_array, ndim, nb, mdist, met, cls_size, folder):
    """runs umap for the specified parameters and returns an array of all euclidean distance pairs between datapoints
    in the umap embeddings """

    UMAPspace = UMAP(n_components=ndim, n_neighbors=nb, min_dist=mdist, metric=met,
                     random_state=0).fit_transform(FP_array)

    UMAPspace = np.nan_to_num(UMAPspace)

    z = linkage(UMAPspace, method="ward")
    cls = fcluster(z, cls_size, criterion="maxclust")
    silhuette_score = silhouette_score(UMAPspace, cls)
    print(nb, mdist, cls_size, silhuette_score)

    title = f"Dim={ndim}_nb={nb}_m-dist={mdist}_cls_size={cls_size}_silscore={silhuette_score}"
    save_umap_fig(UMAPspace, title, cls, folder)
    datasize = FP_array.shape[0]
    norm_reduced_dist = compute_reduced_distances(UMAPspace, datasize)
    return norm_reduced_dist


def true_jaccard_dist(fp_list):
    """returns an array of jaccard distance among all pairs of fps
    Input: list of fps
    Output: np array of jaccard distance pairs among the fps"""

    datasize = len(fp_list)
    sim_array = np.zeros((datasize, datasize))
    for i in range(datasize):
        fp = fp_list[i]

        for j in range(datasize):
            if j != i:
                fp2 = fp_list[j]
                similarity = DataStructs.FingerprintSimilarity(fp, fp2)
                sim_array[i, j] = 1 - similarity

    return sim_array


def save_umap_fig(UMAPspace, title, cls, folder):
    """plots the umap embeddings and saves it in the desired location"""

    fig, ax = plt.subplots()
    plt.scatter(UMAPspace[:, 0], UMAPspace[:, 1], s=8, c=cls, alpha=0.5)
    ax.set_xlabel(r'UMAP1', fontsize=15)
    ax.set_ylabel(r'UMAP2', fontsize=15)
    ax.set_title(title)
    fig.tight_layout()
    filename = "umap_" + title + ".png"
    plt_location = folder.joinpath(filename)
    plt.savefig(plt_location, dpi=120, transparent=True, bbox_inches="tight")
    plt.clf()
    plt.close('all')


def compute_reduced_distances(UMAPspace, datasize):
    """computes the euclidean distance between all pairs of datapoints in the umap embedding"""

    umap_dist_matrix = np.zeros((datasize, datasize))
    for i in range(datasize):
        molecule = UMAPspace[i, :]

        for j in range(datasize):
            if j != i:
                molecule2 = UMAPspace[j, :]
                similarity = distance.euclidean(molecule, molecule2)
                umap_dist_matrix[i, j] = similarity

    reduced_dist_vec = umap_dist_matrix[np.triu_indices(datasize)]
    norm_reduced_dist = reduced_dist_vec / np.max(reduced_dist_vec)
    return norm_reduced_dist


def save_boxplot(original_dist, reduced_dist, title, folder, intervals):
    """plotting boxplot: original distance vs embedded distance"""

    range_list = np.zeros(original_dist.size)
    for i in range(intervals):
        lower_lim = i / intervals
        upper_lim = lower_lim + (1 / intervals)
        if i < (intervals - 1):
            index = np.where(np.logical_and(original_dist >= lower_lim, original_dist < upper_lim))
        else:
            index = np.where(np.logical_and(original_dist >= lower_lim, original_dist <= upper_lim))
        range_list[index] = i

    boxp_df = pd.DataFrame({'original': original_dist.tolist(), 'reduced': reduced_dist.tolist(),
                            'range_list': range_list.tolist()})
    sns.set_theme(style="whitegrid")
    ax = sns.boxplot(x='range_list', y='reduced', data=boxp_df, fliersize=1)
    ax.set_xlabel(r'original distance', fontsize=22)
    ax.set_ylabel(r'embedded distance', fontsize=22)
    filename = "b_" + title + ".png"
    plt_location = folder.joinpath(filename)
    plt.savefig(plt_location, transparent=True, bbox_inches="tight")
    plt.clf()
    plt.close('all')
    print('asd')
