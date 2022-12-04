from copy import deepcopy
from typing import Any, Dict, List, Tuple, Union
import numpy as np

from ase import Atoms

import utils as ut


def rmsd_of_distances(mol1: Atoms, mol2: Atoms) -> float:
  d1 = mol1.get_all_distances()
  d2 = mol2.get_all_distances()
  return np.sqrt(np.mean((d1 - d2)**2))


def rmsd_of_positions(mol1: Atoms, mol2: Atoms) -> float:
  return np.sqrt(np.mean((mol1.positions - mol2.positions)**2))

