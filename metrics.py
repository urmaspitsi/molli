from copy import deepcopy
from typing import Any, Dict, List, Tuple, Union
import numpy as np

from ase import Atoms


def mad_of_positions(mol1: Atoms, mol2: Atoms) -> float:
  '''
    Mean Absolute Distance between position coordinates.
  '''
  return np.mean(np.linalg.norm(mol1.positions - mol2.positions, axis=1))
  #return np.mean(np.abs(mol1.positions - mol2.positions))


def mad_of_distances(mol1: Atoms, mol2: Atoms) -> float:
  '''
    Mean Absolute Difference across all inter-atomic distances.
  '''
  d1 = mol1.get_all_distances()
  d2 = mol2.get_all_distances()
  return np.mean(np.abs(d1 - d2))


def rmsd_of_distances(mol1: Atoms, mol2: Atoms) -> float:
  d1 = mol1.get_all_distances()
  d2 = mol2.get_all_distances()
  return np.sqrt(np.mean((d1 - d2)**2))


def rmsd_of_positions(mol1: Atoms, mol2: Atoms) -> float:
  '''
    Root Mean Squared Distance:
    Square root of the Average of Squared Distances between position coordinates.
  '''
  return np.sqrt(np.mean(np.linalg.norm(mol1.positions - mol2.positions, axis=1)**2))
  #return np.sqrt(np.mean((mol1.positions - mol2.positions)**2))

