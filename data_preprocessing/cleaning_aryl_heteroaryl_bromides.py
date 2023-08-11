### load dataset exported from Reaxys (c)
df = pd.read_csv(path)

### renaming reaxys export column names and changing datatypes for clarity
rename_dict = {
    'Substance Identification: Reaxys Registry Number': 'reaxys_id',
    'Links to Reaxys': 'reaxys_link',
    'Data Count': 'subset_download_id',
    'CAS Registry Number': 'cas_number',
    'Molecular Formula': 'molecular_formula',
    'Molecular Weight': 'molecular_weight',
    'InChI Key': 'inchi_key',
    'Number of Reactions': 'number_reactions',
    'Number of References': 'number_references'
}
datatype_dict = {
    'reaxys_id': int,
    'number_reactions': int,
    'number_references': int
}

df = df.rename(columns=rename_dict)
df = df.fillna(0).astype(datatype_dict)

### filter settings:
max_mw = 750
min_rxn_references = 25

### create canonical smiles
df['mol_objects'] = [MolFromSmiles(str(smiles)) for smiles in df['SMILES'].tolist()]
df['canonical_smiles'] = [(MolToSmiles(mol, canonical=True) if mol is not None else 'none') for mol in df['mol_objects'].tolist()]
df['n_ar_bromine_atoms'] = [(len(mol.GetSubstructMatches(MolFromSmarts('aBr'))) if mol is not None else 0) for mol in df['mol_objects'].tolist()]
df.drop(columns=['mol_objects'], inplace=True)

### remove entries which are not RDKIT-parsable
len_before = len(df)
df = df[~(df['canonical_smiles'] == "none")]
print(f'Removed {len_before - len(df)} entries that were not RDKit-parsable.')

### remove duplicates
len_before = len(df)
df.drop_duplicates('canonical_smiles', inplace=True)
print(f'Removed {len_before - len(df)} duplicates.')

### remove entries with more than one bromine atom
len_before = len(df)
df = df[df['n_ar_bromine_atoms'] == 1]
print(f'Removed {len_before - len(df)} entries that contained more than one bromine atom.')

### remove compounds exceeding the maximum weight
len_before = len(df)
df = df[round(df['molecular_weight']) <= max_mw]
print(f'Removed {len_before - len(df)} entries that were heavier than 750 Da.')

### remove compounds that were used in less than N reactions
len_before = len(df)
df = df[round(df['number_reactions']) >= min_rxn_references]
print(f'Removed {len_before - len(df)} entries that were used in less than {min_rxn_references} reactions.')

### remove compounds that were mentioned in less then N references
len_before = len(df)
df = df[round(df['number_references']) >= min_rxn_references]
print(f'Removed {len_before - len(df)} entries that were referenced in less than {min_rxn_references} papers.')

### remove compounds that contain metals or salts and isotope-marked compounds
len_before = len(df)
removal_characters = ['\[Cd\]', '\[In\]', '\[Ce\]', '\[Gd\]', '\[Os', '\[U\]', '\[K\]', '\[Rb\]', '\[Cs\]', '\[Re\]',
                        '\[Be\]', '\[Hf\]', '\[Bi\]', '\[Ru\]', '\[Li\]', '\[Tl\]', '\[Mg\]', '\[Al\]', '\[Na\]',
                        '\[Pt\]', '\[2H]', '\[3H]', '\[13C]',
                        '\[Mn\]', '\[Ag\]', '\[Sb\]', '\[W\]', '\[Pb\]', '\[Mo', '\[Cu\]', '\[Zr\]', '\[Pd\]', '\[Hg',
                        '\[Te\]', '\[Ti\]', '\[Zn', '\[V\]', '\[Rh\]', '\[Cr\]', '\[Ni\]', '\[Ge', '\[As', '\[Sn',
                        '\[Fe', '\[Co\]', '\[Se', '\[Co', '\[te\]', '\[se\]', '\[Th', '\[Sc', '\[Nb', '\.']
for element in removal_characters:
    df = df[~df['SMILES'].str.contains(element)]
print(f'Removed {len_before - len(df)} compounds containing salts, metals or unnatural isotopes.')


### write cleaned dataset to csv
df.to_csv('cleaned_dataset.csv', index=False)
df[['SMILES','reaxys_id','canonical_smiles', 'cas_number']].to_csv('cleaned_dataset_minimal.csv', index=False)