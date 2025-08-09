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

## Add your files

- [ ] [Create](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#create-a-file) or [upload](https://docs.gitlab.com/ee/user/project/repository/web_editor.html#upload-a-file) files
- [ ] [Add files using the command line](https://docs.gitlab.com/ee/gitlab-basics/add-file.html#add-a-file-using-the-command-line) or push an existing Git repository with the following command:

```
cd existing_repo
git remote add origin https://zivgitlab.uni-muenster.de/ag-glorius/published-paper/standardizing-substrate-selection.git
git branch -M main
git push -uf origin main
```

## Integrate with your tools

- [ ] [Set up project integrations](https://zivgitlab.uni-muenster.de/ag-glorius/published-paper/standardizing-substrate-selection/-/settings/integrations)

## Collaborate with your team

- [ ] [Invite team members and collaborators](https://docs.gitlab.com/ee/user/project/members/)
- [ ] [Create a new merge request](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html)
- [ ] [Automatically close issues from merge requests](https://docs.gitlab.com/ee/user/project/issues/managing_issues.html#closing-issues-automatically)
- [ ] [Enable merge request approvals](https://docs.gitlab.com/ee/user/project/merge_requests/approvals/)
- [ ] [Set auto-merge](https://docs.gitlab.com/ee/user/project/merge_requests/merge_when_pipeline_succeeds.html)

## Test and Deploy

Use the built-in continuous integration in GitLab.

- [ ] [Get started with GitLab CI/CD](https://docs.gitlab.com/ee/ci/quick_start/index.html)
- [ ] [Analyze your code for known vulnerabilities with Static Application Security Testing(SAST)](https://docs.gitlab.com/ee/user/application_security/sast/)
- [ ] [Deploy to Kubernetes, Amazon EC2, or Amazon ECS using Auto Deploy](https://docs.gitlab.com/ee/topics/autodevops/requirements.html)
- [ ] [Use pull-based deployments for improved Kubernetes management](https://docs.gitlab.com/ee/user/clusters/agent/)
- [ ] [Set up protected environments](https://docs.gitlab.com/ee/ci/environments/protected_environments.html)

***

# Editing this README

When you're ready to make this README your own, just edit this file and use the handy template below (or feel free to structure it however you want - this is just a starting point!). Thank you to [makeareadme.com](https://www.makeareadme.com/) for this template.

## Suggestions for a good README
Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.

## Name
Choose a self-explaining name for your project.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
