from copy import deepcopy
from typing import Any, Dict, List, Tuple, Union
import numpy as np
from scipy.spatial import distance_matrix

from ase import Atoms


def average_distance(atoms_obj: Atoms, atom_idxs: List[int]) -> float:
  '''
    Average across all inter-atomic distances.
  '''
  pos = atoms_obj.positions \
    if atom_idxs == None or len(atom_idxs) < 2 \
    else atoms_obj.positions[atom_idxs]
  
  return np.mean(distance_matrix(pos, pos))


def mad_of_positions(mol1: Atoms, mol2: Atoms) -> float:
  '''
    Mean Absolute Distance between position coordinates.
  '''
  return np.mean(np.linalg.norm(mol1.positions - mol2.positions, axis=1))


def mad_of_distances(mol1: Atoms, mol2: Atoms) -> float:
  '''
    Mean Absolute Difference across all inter-atomic distances.
  '''
  d1 = mol1.get_all_distances()
  d2 = mol2.get_all_distances()
  return np.mean(np.abs(d1 - d2))


def mean_abs_coord_diff(mol1: Atoms, mol2: Atoms) -> float:
  '''
    Mean Absolute Difference between position coordinates.
    NB! This is not euclidean distance!
    For euclidean distance use mad_of_positions()
    i.e. mean-absolute-distance
  '''
  return np.mean(np.abs(mol1.positions - mol2.positions))


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


def volume_cube(atoms_obj: Atoms, atom_idxs: List[int]) -> float:
  '''
    Volume of a cube that fits given molecule(atoms) inside.
  '''

  pos = atoms_obj.positions \
    if atom_idxs == None or len(atom_idxs) < 2 \
    else atoms_obj.positions[atom_idxs]
    
  return np.product(np.max(pos, axis=0) - np.min(pos, axis=0))
