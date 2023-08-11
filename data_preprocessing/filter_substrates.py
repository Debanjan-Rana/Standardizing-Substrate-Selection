import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdqueries
from rdkit.Chem.Descriptors import MolWt as MW


class Clean_Smiles():
    """Class for cleaning and filtering substrates"""

    cannonical_smiles: pd.DataFrame
    patt_molwt_processed: pd.DataFrame
    cleaned_smiles_df: pd.DataFrame

    def __init__(self, data_path):
        """ initializes by inputting the folder path of the file containing list of substrate molecules in csv format"""

        self.original_df = pd.read_csv(str(data_path))
        print(f"Original dataset size:{self.original_df.shape[0]}")
        self.cannonical_smiles_df = self.__check_parsing_error()

    def __check_parsing_error(self) -> pd.DataFrame:
        """checking if any molecules/smiles have rdkit parsing error"""

        cannonical_smiles_list = []
        error_parsing_list = []
        drop_index_list = []
        for i in range(self.original_df.shape[0]):
            smiles = self.original_df['SMILES'][i]
            try:
                cannonical_smiles_list.append(
                    Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=False, canonical=True))
            except:
                print('error_processing_molecule')
                error_parsing_list.append(smiles)
                drop_index_list.append(i)
        print(f"Molecules without rdkit parsing error: {len(cannonical_smiles_list)}")

        self.error_parsing = error_parsing_list
        cannonical_smiles_df = pd.DataFrame({'SMILES': cannonical_smiles_list})
        return cannonical_smiles_df

    def clean_carbon_molwt(self, smiles_df, n_Carbon, l_mol_wt, u_mol_wt) -> pd.DataFrame:
        """Cleaning for molecules which has less than 2 Carbon atom and molecular weight criteria"""

        q = rdqueries.AtomNumEqualsQueryAtom(6)
        match_molecules = []
        excluded_smiles = []
        drop_index_list = []
        for i in range(smiles_df.shape[0]):
            smiles = smiles_df['SMILES'][i]
            mol = Chem.MolFromSmiles(smiles)

            if len(mol.GetAtomsMatchingQuery(q)) > n_Carbon and u_mol_wt > MW(mol) > l_mol_wt:
                match_molecules.append(Chem.MolToSmiles(mol))
                # print(smiles)
            else:
                excluded_smiles.append(Chem.MolToSmiles(mol))
                drop_index_list.append(i)

        self.excluded_smiles = excluded_smiles
        print("")
        print(f"Molecules having atleast 2 Carbon atoms and weight less than {u_mol_wt} and greater than {l_mol_wt}: {len(match_molecules)}")
        matched_pattern_df = pd.DataFrame({'SMILES': match_molecules})
        return matched_pattern_df

    def cleaning_salts(self, smiles_df, pattern) -> pd.DataFrame:
        """Cleaning for salt and considering counterparts that contain the desired substructure. The substructure
        needs to be specified as SMARTS. e.g., for alkenes [#6]=[#6]"""

        alkene_fragemnt = []
        non_salts = []
        pattern_mol = Chem.MolFromSmarts(pattern)
        count = 0
        for i in range(smiles_df.shape[0]):
            smiles = smiles_df['SMILES'][i]
            if '.' in smiles:
                count += 1
                ions = smiles.split('.')
                for ion in ions:
                    ion_mol = Chem.MolFromSmiles(ion)
                    if ion_mol.HasSubstructMatch(pattern_mol):
                        alkene_fragemnt.append(Chem.MolToSmiles(ion_mol))

            else:
                non_salts.append(smiles)
        cleaned_smiles_list = non_salts + alkene_fragemnt
        cleaned_smiles_df = pd.DataFrame({'SMILES': cleaned_smiles_list})
        print("")
        print(f"Total molecules after cleaning salts: {cleaned_smiles_df.shape[0]}")
        return cleaned_smiles_df

    def screen_pattern(self, smiles_df, pattern, type='include') -> pd.DataFrame:
        """Excluding incompatible functional groups/ substructures or including compatible substructures depending on
        the reactivity knowledge. """

        matched_pattern = []
        for smiles in smiles_df['SMILES']:
            mol = Chem.MolFromSmiles(smiles)
            if type == 'include':
                pattern_mol = Chem.MolFromSmarts(pattern)
                if mol.HasSubstructMatch(pattern_mol):
                    matched_pattern.append(Chem.MolToSmiles(mol))
                else:
                    pass
                    # print(smiles)

            """includes molecules that matches either of the specified smarts substructures"""
            if type == 'multi':
                pattern_mol1 = Chem.MolFromSmarts(pattern[0])
                pattern_mol2 = Chem.MolFromSmarts(pattern[1])
                pattern_mol3 = Chem.MolFromSmarts(pattern[2])
                pattern_mol4 = Chem.MolFromSmarts(pattern[3])
                if mol.HasSubstructMatch(pattern_mol1) or mol.HasSubstructMatch(pattern_mol2) or mol.HasSubstructMatch(
                        pattern_mol3) or mol.HasSubstructMatch(pattern_mol4):
                    matched_pattern.append(Chem.MolToSmiles(mol))
                else:
                    pass
                    # print(smiles)

            """excludes from the dataframe doing the opposite of include operation"""
            if type == 'exclude':
                pattern_mol = Chem.MolFromSmarts(pattern)
                if mol.HasSubstructMatch(pattern_mol):
                    pass
                    # print(smiles)
                else:
                    matched_pattern.append(Chem.MolToSmiles(mol))

        if type == 'exclude':
            print("")
            print(f"Molecules without {pattern}: {len(matched_pattern)}")
        else:
            print("")
            print(f"Molecules having {pattern}: {len(matched_pattern)}")

        matched_pattern_df = pd.DataFrame({'SMILES': matched_pattern})
        return matched_pattern_df

    def remove_duplicates(self, smiles_df) -> pd.DataFrame:
        """removing duplicates in the filtered list"""

        clean_df = smiles_df.drop_duplicates(subset='SMILES')
        clean_df.reset_index(drop=True, inplace=True)
        print(f"Final_dataset_after_removing_duplicates: {clean_df.shape[0]}")
        return clean_df


RDLogger.DisableLog('rdApp.*')
# folder_path = Path.cwd()

path = 'raw_data/exported_alkenes.csv'
smiles_preprocessor = Clean_Smiles(path)
alkenes_filtered_df = smiles_preprocessor.screen_pattern(smiles_preprocessor.cannonical_smiles_df, '[#6]=[#6]', type='include')
salts_filtered_df = smiles_preprocessor.cleaning_salts(alkenes_filtered_df, '[#6]=[#6]')
duplicates_filtered_df = smiles_preprocessor.remove_duplicates(salts_filtered_df)

mw_wt_filtered_df = smiles_preprocessor.clean_carbon_molwt(duplicates_filtered_df, n_Carbon=1, l_mol_wt=0, u_mol_wt=700)

fg_filtered_df1 = smiles_preprocessor.screen_pattern(mw_wt_filtered_df, '[NH2]~*', type='exclude')
patterns = ['[CH2]=[CH0]', '[CH2]=[CH1]', '[CH1]=[CH1]', '[CH0]=[CH1]']
final_alkenes = smiles_preprocessor.screen_pattern(fg_filtered_df1, patterns, type='multi')
final_alkenes.to_csv('raw_data/amine_tetra_cleaned.csv', index=False)

print('wego')

