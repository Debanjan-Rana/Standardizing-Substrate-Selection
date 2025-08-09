# Standardizing Substrate Selection

With over 10,000 new reaction protocols arising every year, only a handful of these procedures transition from academia to application. A major reason for this gap stems from the lack of comprehensive knowledge about a reaction’s scope — i.e., to which substrates the protocol can or cannot be applied.

Even though chemists invest substantial effort to assess the scope of new protocols, the resulting scope tables involve significant biases, reducing their expressiveness.

Herein we report a **standardized substrate selection strategy** designed to mitigate these biases and evaluate the applicability, as well as the limits, of any chemical reaction.

- **Unsupervised learning** is utilized to map the chemical space of industrially relevant molecules.
- Potential substrate candidates are projected onto this universal map, enabling the selection of a structurally diverse set of substrates with optimal relevance and coverage.
- By testing our methodology on different chemical reactions, we demonstrate its effectiveness in finding general reactivity trends using only a few highly representative examples.

The developed methodology empowers chemists to showcase the **unbiased applicability** of novel methodologies, facilitating their practical applications.  
We hope this work triggers interdisciplinary discussions about biases in synthetic chemistry, leading to improved data quality.

**📄 Link to the paper:** [ACS Central Science](https://pubs.acs.org/doi/10.1021/acscentsci.3c01638)
**📄 Link to the web-interface:** [Pharma Scope](https://pharmascope.uni-muenster.de/)
---

## Implementation Guidelines

### Files Overview
- **`standardized_scope.txt`** — All dependencies required to run the Python scripts.
- **`scope_generator.py`** — Main workflow execution file.
- **`substrate_selection.py`** — Contains the three main classes.
- **`settings.json`** — All workflow parameters and settings.

### Class Descriptions
#### 1. `Fingerprints`
- Generates and stores specified molecular fingerprints and SMILES for all molecules.

#### 2. `Chemspace_mapper`
- Runs UMAP dimensionality reduction.
- Performs hierarchical agglomerative clustering.
- Selects substrates based on the specified selection strategy.

#### 3. `Scope_manager`
- Wrapper class that incorporates the functionalities of `Fingerprints` and `Chemspace_mapper`.
- Reads the settings from `settings.json`.
- Loads datasets and saves obtained results.

---

## `settings.json` Parameters

```json
{
  "dataset": ["path/to/drug_dataset_folder"],
  "test_dataset": ["path/to/substrate_dataset_folder"],
  "fp_settings": {
    "type": "ECFP", 
    "radius": 2,
    "nBits": 1024
  },
  "umap_settings": {
    "n_neighbors": 15,
    "min_dist": 0.1
  },
  "n_clusters": 10,
  "additional_settings": {
    "load_model": false,
    "model_path": "model_name.sav",
    "exp_name": "experiment_folder_name",
    "topn_mol": 3,
    "selection_strategy": "centre-based"
  }
}
```

## Parameter Details

- **dataset** *(list)*: Folder location for the drug dataset.
- **test_dataset** *(list)*: Folder location for the respective substrate dataset.
- **fp_settings** *(dict)*: Specifies the fingerprint type.  
  - Supported: `ECFP`, `MACCS`.  
  - For `ECFP`, also specify the `radius` and number of bits (`nBits`).
- **umap_settings** *(dict)*: UMAP parameter configuration.
- **n_clusters** *(int)*: Number of clusters for hierarchical agglomerative clustering.

**additional_settings** *(dict)*:
- **load_model** *(bool)*:  
  - `true`: Use a pre-trained model.  
  - `false`: Train a new model with the given UMAP parameters.
- **model_path** *(str)*: Path to the model file (required if `load_model` is `true`, e.g., `model_name.sav`).
- **exp_name** *(str)*: Folder name for saving experiment results.
- **topn_mol** *(int)*: Number of substrates to select from each cluster, based on the selection strategy.
- **selection_strategy** *(str)*:  
  - `centre-based`: Selects the substrate closest to the cluster centroid.  
  - `similarity-based`: Selects the substrate dtructurally most similar to all drugs in that cluster.
