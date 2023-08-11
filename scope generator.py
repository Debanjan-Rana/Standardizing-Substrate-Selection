from rdkit import RDLogger
from substrate_selection import Scope_Manager

RDLogger.DisableLog('rdApp.*')
space = Scope_Manager()
space.save_all_results()
