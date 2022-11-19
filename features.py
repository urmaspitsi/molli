from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple, Union

from ase import Atoms


@dataclass
class Feature():
  label: str

  def calculate_value(self, atoms_obj: Atoms) -> float:
    raise NotImplementedError

  def get_info(self, atoms_obj: Atoms) -> str:
    raise NotImplementedError


@dataclass
class Bond(Feature):
  atom_idx1: int
  atom_idx2: int

  def calculate_value(self, atoms_obj: Atoms) -> float:
    return atoms_obj.get_distance(self.atom_idx1, self.atom_idx2)

  def get_info(self, atoms_obj: Atoms) -> str:
    symbols = atoms_obj.get_chemical_symbols()
    n1 = f"{symbols[self.atom_idx1]}({self.atom_idx1})"
    n2 = f"{symbols[self.atom_idx2]}({self.atom_idx2})"
    return f"Distance between atoms {n1}-{n2}"


@dataclass
class Angle(Feature):
  atom_idx1: int
  atom_idx2: int
  atom_idx3: int

  def calculate_value(self, atoms_obj: Atoms) -> float:
    return atoms_obj.get_angle(self.atom_idx1, self.atom_idx2, self.atom_idx3)

  def get_info(self, atoms_obj: Atoms) -> str:
    symbols = atoms_obj.get_chemical_symbols()
    n1 = f"{symbols[self.atom_idx1]}({self.atom_idx1})"
    n2 = f"{symbols[self.atom_idx2]}({self.atom_idx2})"
    n3 = f"{symbols[self.atom_idx3]}({self.atom_idx3})"
    return f"Bending angle between atoms {n1}-{n2}-{n3}"


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
    n1 = f"{symbols[self.atom_idx1]}({self.atom_idx1})"
    n2 = f"{symbols[self.atom_idx2]}({self.atom_idx2})"
    n3 = f"{symbols[self.atom_idx3]}({self.atom_idx3})"
    n4 = f"{symbols[self.atom_idx4]}({self.atom_idx4})"
    return f"Torsion angle between atoms {n1}-{n2}-{n3}-{n4}"


@dataclass
class CellVolume(Feature):
  atom_idxs: List[int] = field(default_factory=list)

  def calculate_value(self, atoms_obj: Atoms) -> float:
    return atoms_obj.get_volume()
