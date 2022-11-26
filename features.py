from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple, Union
import numpy as np
from scipy.spatial import distance_matrix

from ase import Atoms


@dataclass
class Feature():
  label: str

  def calculate_value(self, atoms_obj: Atoms) -> float:
    raise NotImplementedError

  def get_info(self, atoms_obj: Atoms) -> str:
    raise NotImplementedError

  def get_units(self) -> str:
    raise NotImplementedError

@dataclass
class Distance(Feature):
  atom_idx1: int
  atom_idx2: int

  def calculate_value(self, atoms_obj: Atoms) -> float:
    return atoms_obj.get_distance(self.atom_idx1, self.atom_idx2)

  def get_info(self, atoms_obj: Atoms) -> str:
    symbols = atoms_obj.get_chemical_symbols()
    n1 = f"{symbols[self.atom_idx1]}({self.atom_idx1 + 1})"
    n2 = f"{symbols[self.atom_idx2]}({self.atom_idx2 + 1})"
    return f"Distance between atoms {n1}-{n2}"

  def get_units(self) -> str:
    return "Angstroms"

@dataclass
class Angle(Feature):
  atom_idx1: int
  atom_idx2: int
  atom_idx3: int

  def calculate_value(self, atoms_obj: Atoms) -> float:
    return atoms_obj.get_angle(self.atom_idx1, self.atom_idx2, self.atom_idx3)

  def get_info(self, atoms_obj: Atoms) -> str:
    symbols = atoms_obj.get_chemical_symbols()
    n1 = f"{symbols[self.atom_idx1]}({self.atom_idx1 + 1})"
    n2 = f"{symbols[self.atom_idx2]}({self.atom_idx2 + 1})"
    n3 = f"{symbols[self.atom_idx3]}({self.atom_idx3 + 1})"
    return f"Bending angle between atoms {n1}-{n2}-{n3}"

  def get_units(self) -> str:
    return "Degrees"

@dataclass
class Dihedral(Feature):
  atom_idx1: int
  atom_idx2: int
  atom_idx3: int
  atom_idx4: int

  def calculate_value(self, atoms_obj: Atoms) -> float:
    return atoms_obj.get_dihedral(self.atom_idx1, self.atom_idx2, self.atom_idx3, self.atom_idx4)

  def get_info(self, atoms_obj: Atoms) -> str:
    symbols = atoms_obj.get_chemical_symbols()
    n1 = f"{symbols[self.atom_idx1]}({self.atom_idx1 + 1})"
    n2 = f"{symbols[self.atom_idx2]}({self.atom_idx2 + 1})"
    n3 = f"{symbols[self.atom_idx3]}({self.atom_idx3 + 1})"
    n4 = f"{symbols[self.atom_idx4]}({self.atom_idx4 + 1})"
    return f"Torsion angle between atoms {n1}-{n2}-{n3}-{n4}"

  def get_units(self) -> str:
    return "Degrees"

@dataclass
class Volume(Feature):
  atom_idxs: List[int] = field(default_factory=list)

  def calculate_value(self, atoms_obj: Atoms) -> float:
    pos = atoms_obj.positions \
      if self.atom_idxs == None or len(self.atom_idxs) < 2 \
      else atoms_obj.positions[self.atom_idxs]
      
    return np.product(np.max(pos, axis=0) - np.min(pos, axis=0))

  def get_info(self, atoms_obj: Atoms) -> str:
    return "Volume"

  def get_units(self) -> str:
    return "Cubic Angstroms"


@dataclass
class AverageDistance(Feature):
  atom_idxs: List[int] = field(default_factory=list)

  def calculate_value(self, atoms_obj: Atoms) -> float:
    pos = atoms_obj.positions \
      if self.atom_idxs == None or len(self.atom_idxs) < 2 \
      else atoms_obj.positions[self.atom_idxs]
    
    return np.mean(distance_matrix(pos, pos))

  def get_info(self, atoms_obj: Atoms) -> str:
    return "Average distance between atoms"

  def get_units(self) -> str:
    return "Angstroms"
