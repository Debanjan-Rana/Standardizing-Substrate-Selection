from rdkit import Chem
from copy import deepcopy
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from typing import List, Dict, Union
from rdkit import RDLogger, DataStructs
from umap import UMAP
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import json

from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score
from scipy.spatial import distance
import pickle

print("Imports completed")

class Fingerprints():
    type: str
    radius: int
    n_bits: int
    raw_fps: list
    smiles: List[str]
    fp_dataframe: pd.DataFrame

    def __init__(self, mol_list: List[Chem.rdchem.Mol], fp_settings: Dict[str, Union[str, int,]]):
        self.type = fp_settings.get("fp_type")
        self.radius = fp_settings.get("fp_radius", 2)
        self.n_bits = fp_settings.get("nr_bits", 1024)
        self.raw_fps = self.__mols_to_fingerprints(mol_list)
        self.smiles = self.__get_smiles(mol_list)
        self.fp_dataframe = self.__construct_fp_dataframe()

    def __mols_to_fingerprints(self, mol_list: List[Chem.Mol]) -> list:
        """generates specified molecular fingerprints from mol objects"""
        if self.type == "MACCS":
            fingerprints = [Chem.MACCSkeys.GenMACCSKeys(mol) for mol in mol_list]
            print(f"MACCSkeys fingerprints generated.")
        elif self.type == "ECFP":
            fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, self.radius, nBits=self.n_bits) for mol in
                            mol_list]
            print(f"ECFP fingerprints with radius {self.radius} and {self.n_bits} bits generated.")
        else:
            print("No fingerprint type given.")
            fingerprints = [Chem.MACCSkeys.GenMACCSKeys(mol) for mol in mol_list]
            print("Default (MACCS keys) used.")
            print("Specify the fingerprint type by the fp_type variable")
        return fingerprints

    def __get_smiles(self, mol_list: List[Chem.rdchem.Mol]) -> List[str]:
        """converts mol objects to smiles"""
        smiles_list = []
        for mol in mol_list:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
            smiles_list.append(smiles)
        return smiles_list

    def __construct_fp_dataframe(self):
        fingerprint_arr = np.array(self.raw_fps)
        fp_dataframe = pd.DataFrame(fingerprint_arr)
        return fp_dataframe


class Chemspace_mapper():
    """Performs umap dimensionality reduction, clustering and substrate selection based on the specified
    settings and input data."""
    drug_embedding: pd.DataFrame
    substrate_embedding: np.array
    n_clusters: int
    silhuette_scores: float
    drug_cluster_ids: np.array
    drug_cluster_mean: np.array
    selected_substrates_df: pd.DataFrame
    simi_selected_substrates_df: pd.DataFrame
    substrate_distance_dict: Dict[int, pd.DataFrame]

    def __init__(self, drug_fingerprints: Fingerprints, subs_fingerprints: Fingerprints,
                 workflow_settings: Dict, additional_settings):

        self.umap_settings = workflow_settings.get("umap_settings")
        self.n_clusters = workflow_settings.get("n_clusters")
        self.fp_settings = workflow_settings.get("fp_settings")
        self.load_model = additional_settings.get("load_model")
        self.model_name = additional_settings.get("model_path")
        self.topn = additional_settings.get("topn_mol")
        self.selection_strategy = additional_settings.get("selection_strategy")
        self.drug_smiles = drug_fingerprints.smiles
        self.substrate_smiles = subs_fingerprints.smiles
        self.drug_embedding = self.__reduce_fingerprints(drug_fingerprints, subs_fingerprints)

    def select_scope(self):
        """clustering drug embeddings"""
        n_clusters = self.n_clusters
        z = linkage(self.drug_embedding, method="ward")
        cls = fcluster(z, n_clusters, criterion="maxclust")
        silhuette_score = silhouette_score(self.drug_embedding, cls)

        self.silhuette_scores = silhuette_score
        self.drug_cluster_ids = deepcopy(cls)
        self.drug_cluster_mean = self.__compute_cluster_mean(n_clusters, cls)

        print("")
        print("selecting substrates for scope")
        self.selected_substrates_df = self.__central_substrates(self.substrate_embedding, self.drug_cluster_mean, n_clusters)
        if self.selection_strategy != "centre-based":
            self.simi_selected_substrates_df = self.__similarity_selection()

    def __reduce_fingerprints(self, drug_fingerprints: Fingerprints, subs_fingerprints: Fingerprints) -> pd.DataFrame:
        """trains umap model on drug fingerprints and then projects substrates with the trained model.
        If load_model is true then a pretrained model is loaded otherwise a new model is trained"""

        dims = self.umap_settings.get("dimensions")
        neighbours = self.umap_settings.get("neighbours")
        min_dist = self.umap_settings.get("m_dist")
        metric = self.umap_settings.get("metric")
        print("")
        print("training UMAP model to map drug fingerprints")
        print(f"Settings: {dims} dims., {neighbours} neighbours, {min_dist} min. distance.")
        drug_fp_dataframe = drug_fingerprints.fp_dataframe

        if self.load_model:
            model_path = Path.cwd().joinpath("trained_models", self.model_name)
            UMAPspace = pickle.load((open(model_path, 'rb')))
        else:
            UMAPspace = UMAP(n_components=dims, n_neighbors=neighbours, min_dist=min_dist,
                             metric=metric, random_state=0).fit(drug_fp_dataframe)
        indexes = drug_fp_dataframe.index
        drug_embedding = pd.DataFrame(UMAPspace.embedding_, index=indexes)
        print('UMAP_done')

        print("")
        print('projecting substrates with trained drug model')
        subs_fp_dataframe = subs_fingerprints.fp_dataframe
        self.substrate_embedding = UMAPspace.transform(subs_fp_dataframe)

        return drug_embedding

    def __similarity_selection(self) -> pd.DataFrame:
        hits = []
        alk_cluster_idx = []
        drug_df = pd.DataFrame({'drugs': self.drug_smiles})
        alk_df = pd.DataFrame({'alk': self.substrate_smiles})

        for i in range(self.n_clusters):
            sim_list = []
            subs_in_cls = alk_df.loc[self.substrate_distance_dict[i]['molecules'].to_list()]
            assgn_alk_mol_list = [Chem.MolFromSmiles(smiles) for smiles in subs_in_cls['alk']]
            assgn_alk_fp_list = [AllChem.GetMorganFingerprintAsBitVect(mol, self.fp_settings['fp_radius'], self.fp_settings['nr_bits']) for mol in assgn_alk_mol_list]

            curr_cls = i + 1
            drugs_in_cls = np.where(self.drug_cluster_ids == curr_cls)[0]
            drugs_in_cls_list = drugs_in_cls.tolist()
            assgn_drug = drug_df.loc[drugs_in_cls_list]
            assgn_drug_mol_list = [Chem.MolFromSmiles(smiles) for smiles in assgn_drug['drugs']]
            assgn_drug_fp_list = [AllChem.GetMorganFingerprintAsBitVect(mol, self.fp_settings['fp_radius'], self.fp_settings['nr_bits']) for mol in assgn_drug_mol_list]

            for alk_fp in assgn_alk_fp_list:
                sm_store = 0
                for drug_fp in assgn_drug_fp_list:
                    sm_morgan = DataStructs.FingerprintSimilarity(alk_fp, drug_fp)
                    sm_store = sm_store + sm_morgan
                avg_sim_per_fp = sm_store / len(assgn_drug_fp_list)
                sim_list.append(avg_sim_per_fp)
            sim_arr = np.array(sim_list)

            num_topn_sim_alk = self.topn
            minimum_mols = len(sim_arr)

            if minimum_mols > num_topn_sim_alk:
                sorted_arr = np.argsort(sim_arr)[::-1][0:num_topn_sim_alk]
                for sim_alk in range(num_topn_sim_alk):
                    hits.append(subs_in_cls['alk'].to_list()[sorted_arr[sim_alk]])
                    alk_cluster_idx.append(i)
            else:
                sorted_arr = np.argsort(sim_arr)[::-1][0:minimum_mols]
                for sim_alk in range(minimum_mols):
                    hits.append(subs_in_cls['alk'].to_list()[sorted_arr[sim_alk]])
                    alk_cluster_idx.append(i)

        simi_selected_subs_df = pd.DataFrame({'smiles': hits, 'substrate_cls_id': alk_cluster_idx})
        return simi_selected_subs_df

    def __compute_cluster_mean(self, n_clusters, cls) -> np.array:
        """returns an array with the mean of each drug cluster"""

        cluster_mean = np.zeros((n_clusters, self.drug_embedding.shape[1]))
        drug_cluster_ids = cls
        grps = self.drug_embedding.groupby(drug_cluster_ids)
        for single_cluster in range(n_clusters):
            cls_mol_indices = grps.indices[single_cluster + 1]
            # computing mean for each cluster
            single_cluster_mean = self.drug_embedding.iloc[cls_mol_indices].mean(0)
            cluster_mean[single_cluster, :] = single_cluster_mean

        return cluster_mean

    def __central_substrates(self, substrate_embedding, drug_cluster_mean, n_clusters) -> pd.DataFrame:
        """selects substrates lying closest to the center of each drug cluster"""
        substrate_cls_ids, substrates_distance_list = self.__assign_substrate_cls(drug_cluster_mean,
                                                                                  substrate_embedding, n_clusters)
        center_substrates_df = self.__compute_central_substrates(substrate_cls_ids, substrates_distance_list, n_clusters)
        return center_substrates_df

    def __assign_substrate_cls(self, drug_cluster_mean, substrate_embedding, n_clusters):
        """assigns substrates to drug clusters based on the proximity to the respective cluster centers and returns a
        list of substrate_cluster_ids and a list of substrate distances to the cluster centers"""
        subs_cls_list = []
        subs_distance_list = []
        for i in range(len(self.substrate_smiles)):
            subs_embedding = substrate_embedding[i]
            subs_dist_all_cluster = [distance.euclidean(subs_embedding, drug_cluster_mean[j]) for j in
                                     range(n_clusters)]  # computing distance from all drug cluster centers
            subs_cluster_num = subs_dist_all_cluster.index(min(subs_dist_all_cluster))
            subs_distance = min(subs_dist_all_cluster)
            subs_cls_list.append(subs_cluster_num)
            subs_distance_list.append(subs_distance)
        return subs_cls_list, subs_distance_list

    def __compute_central_substrates(self, substrate_cls_ids, substrates_distance_list, n_clusters) -> pd.DataFrame:
        """returns dataframe with topn molecules closest to the cluster centers. {smiles, mol_ids, cluster_ids}"""
        topn_mol_cls_id = []
        topn_substrates = []
        substrate_distance_dict = {}

        subs_cluster_arr = np.array(substrate_cls_ids)
        subs_distance_arr = np.array(substrates_distance_list)

        for k in range(n_clusters):
            substrate_cls = np.where(subs_cluster_arr == k)[0]
            substrate_dist = subs_distance_arr[substrate_cls]
            single_cls_subs_df = pd.DataFrame({'molecules': substrate_cls, 'mol_dist_from_center': substrate_dist})
            substrate_distance_dict.update({k: single_cls_subs_df})
            topn_subs_single_cluster = self.__topn_mol(single_cls_subs_df)
            topn_substrates = topn_substrates + topn_subs_single_cluster

            mol_cls_ids = [k]*len(topn_subs_single_cluster)
            topn_mol_cls_id = topn_mol_cls_id + mol_cls_ids

        # dict of dataframes containing substrate_ids and distance to cluster center
        self.substrate_distance_dict = substrate_distance_dict

        selected_substrate_smiles = self.__index_to_smiles(topn_substrates, self.substrate_smiles)
        selected_substrates_df = pd.DataFrame({'smiles': selected_substrate_smiles, 'substrate_id': topn_substrates, 'substrate_cls_id': topn_mol_cls_id})
        return selected_substrates_df

    def __topn_mol(self, cls_dist_dataframe: pd.DataFrame) -> List[int]:
        """sorts substrates based on cluster center distance and returns a list of topn closest central molecules"""
        sorted_cls_dist_dataframe = cls_dist_dataframe.sort_values(by='mol_dist_from_center')
        top_cluster_molecules = sorted_cls_dist_dataframe.iloc[0:self.topn]['molecules'].to_list()
        return top_cluster_molecules

    def __index_to_smiles(self, center_molecules_list, smiles_list):
        """returns a list of smiles for the corresponding substrate indices"""
        center_smiles = [smiles_list[index] for index in center_molecules_list]
        return center_smiles


class Scope_Manager():
    settings: Dict[str, Union[str, int, list, dict]]
    dataset: List[Chem.rdchem.Mol]
    drug_fingerprints: Fingerprints
    subs_fingerprints: Fingerprints
    drug_substrate_map: Chemspace_mapper
    test_data: List[Chem.rdchem.Mol]

    def __init__(self):
        self.settings = self.__read_settings()
        data_settings = self.settings.get("dataset")
        self.dataset = self.__read_dataset(data_settings)
        fp_settings = self.settings.get("fp_settings")
        self.drug_fingerprints = Fingerprints(self.dataset, fp_settings)
        self.additional_settings = self.settings.get("additional_settings")
        self.test_data = self.__read_test_dataset()
        self.subs_fingerprints = Fingerprints(self.test_data, fp_settings)
        self.drug_substrate_map = self.__build_drug_substrate_map()

    def __build_drug_substrate_map(self):
        workflow_settings = self.settings
        drug_fingerprints = self.drug_fingerprints
        subs_fingerprints = self.subs_fingerprints
        addn_settings = self.settings.get("additional_settings")
        chemspace_map = Chemspace_mapper(drug_fingerprints, subs_fingerprints, workflow_settings, addn_settings)
        chemspace_map.select_scope()
        return chemspace_map

    # Settings reading
    def __read_settings(self) -> Dict[str, Union[str, int, list]]:
        settings_path = Path.cwd().joinpath("settings", "settings.json")
        with open(settings_path, "r") as settings_file:
            setting_dict = json.load(settings_file)
        return setting_dict

    def __read_dataset(self, path_list) -> List[Chem.rdchem.Mol]:
        data_paths = self.__construct_data_paths(path_list)
        print("")
        print("Reading Dataset!")
        print("Reader implemented for the following file-types: .sdf, .csv")
        print("")
        mol_list = []
        for data_path in data_paths:
            if data_path.suffix == ".sdf":
                mol_list = self.__read_dataset_from_sdf(data_path)
                print(f"Dataset {data_path.name} read!")
                print(f"{len(mol_list)} molecules found.")

            elif data_path.suffix == ".csv":
                single_mol_list = self.__read_dataset_from_csv(data_path)
                print(f"Dataset {data_path.name} read!")
                mol_list = mol_list + single_mol_list
                print(f"{len(mol_list)} molecules found.")
            else:
                print(f"Datatype {data_path.suffix} not found.")
        print("")
        return mol_list

    def __read_test_dataset(self) -> List[Chem.rdchem.Mol]:
        test_data_settings = self.settings.get("test_dataset")
        test_data = self.__read_dataset(test_data_settings)
        return test_data

    def __construct_data_paths(self, path_list) -> List[Path]:
        folder_path = Path.cwd()
        for element in path_list:
            folder_path = folder_path.joinpath(element)
        data_paths = sorted(folder_path.glob("**/*"))
        return data_paths

    def __read_dataset_from_sdf(self, data_path: Path) -> List[Chem.rdchem.Mol]:
        suppl = Chem.SDMolSupplier(str(data_path))
        mol_list = []
        for i, mol in enumerate(suppl):
            try:
                mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True))
                mol_list.append(mol)
            except:
                pass
        return mol_list

    def __read_dataset_from_csv(self, data_path: Path) -> List[Chem.rdchem.Mol]:
        mol_data = pd.read_csv(str(data_path))['SMILES']
        single_mol_list = []
        try:
            single_mol_list = [Chem.MolFromSmiles(smiles) for smiles in mol_data]
        except:
            pass
        return single_mol_list

    def save_all_results(self):
        """saves the selected substrates and generates plots of the umap embeddings"""
        fp_settings = self.settings.get("fp_settings")
        exp_name = self.additional_settings.get("exp_name")
        subfolder_name = f"{exp_name}_{fp_settings.get('fp_type')}_{fp_settings.get('fp_radius')}_{fp_settings.get('nr_bits')}"
        folder = Path.cwd().joinpath("results", subfolder_name)
        folder.mkdir(exist_ok=True)
        drug_substrate_map = self.drug_substrate_map
        self.__plot_embeddings(drug_substrate_map, folder)
        self.__save_selected_molecules(drug_substrate_map, folder)

    def __save_selected_molecules(self, drug_substrate_map, folder):
        """saves the selected substrates information"""
        final_substrate_df = drug_substrate_map.selected_substrates_df
        substrate_filename = "centre_based.csv"
        file_location = folder.joinpath(substrate_filename)
        final_substrate_df.to_csv(file_location)

        if self.additional_settings.get("selection_strategy") != "centre-based":
            simi_final_substrate_df = drug_substrate_map.simi_selected_substrates_df
            substrate_filename = "similarity_based.csv"
            file_location = folder.joinpath(substrate_filename)
            simi_final_substrate_df.to_csv(file_location)

    def __plot_embeddings(self, drug_substrate_map: Chemspace_mapper, folder: Path):
        """returns scatter plot of the overlaid drug and substrate embeddings"""
        reduced_space_two_d = drug_substrate_map.drug_embedding.loc[:, [0, 1]]
        x = reduced_space_two_d.loc[:, 0]
        y = reduced_space_two_d.loc[:, 1]
        cluster_numbers = drug_substrate_map.drug_cluster_ids
        score = drug_substrate_map.silhuette_scores

        title = f"HAC_nb={drug_substrate_map.umap_settings['neighbours']}_m-dist={drug_substrate_map.umap_settings['m_dist']}_score={score}"
        filename = title + ".png"
        plt_location = folder.joinpath(filename)

        fig, ax = plt.subplots()
        plt.scatter(x=x, y=y, c=cluster_numbers, cmap='rainbow', alpha=0.4)

        x_test = drug_substrate_map.substrate_embedding[:, 0]
        y_test = drug_substrate_map.substrate_embedding[:, 1]
        plt.scatter(x=x_test, y=y_test, s=8, alpha=0.8, marker="D")

        ax.set_xlabel(r'UMAP1', fontsize=15)
        ax.set_ylabel(r'UMAP2', fontsize=15)
        ax.set_title(title)
        fig.tight_layout()
        # plt.show()
        plt.savefig(plt_location, dpi=120, transparent=True, bbox_inches="tight")
        plt.clf()
        plt.close('all')


# RDLogger.DisableLog('rdApp.*')
# space = Scope_Manager()
# space.save_all_results()

# print("Woop")
