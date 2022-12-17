from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple, Union
import numpy as np

from ase import Atoms
import metrics as ms

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
    return f"Distance {n1}-{n2}"

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
    return f"Angle {n1}-{n2}-{n3}"

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
    return f"Torsion angle {n1}-{n2}-{n3}-{n4}"

  def get_units(self) -> str:
    return "Degrees"

@dataclass
class Volume(Feature):
  atom_idxs: List[int] = field(default_factory=list)

  def calculate_value(self, atoms_obj: Atoms) -> float:
    return ms.volume_cube(
                          atoms_obj=atoms_obj,
                          atom_idxs=self.atom_idxs
                          )

  def get_info(self, atoms_obj: Atoms) -> str:
    return "Volume"

  def get_units(self) -> str:
    return "Cubic Angstroms"


@dataclass
class AverageDistance(Feature):
  atom_idxs: List[int] = field(default_factory=list)

  def calculate_value(self, atoms_obj: Atoms) -> float:
    return ms.average_distance(
        atoms_obj=atoms_obj,
        atom_idxs=self.atom_idxs
      )

  def get_info(self, atoms_obj: Atoms) -> str:
    return "Average distance between atoms"

  def get_units(self) -> str:
    return "Angstroms"


def create_angle(
    atom_idxs: Tuple[int],
    mol: Atoms=None) -> Angle:

  idx1, idx2, idx3 = atom_idxs
  res = Angle(
              label=f"{idx1 + 1}-{idx2 + 1}-{idx3 + 1}",
              atom_idx1=idx1,
              atom_idx2=idx2,
              atom_idx3=idx3
            )
  if mol != None:
    res.label = res.get_info(atoms_obj=mol)

  return res


def create_dihedral(
    atom_idxs: Tuple[int],
    mol: Atoms=None) -> Dihedral:

  idx1, idx2, idx3, idx4 = atom_idxs
  res = Dihedral(
              label=f"{idx1 + 1}-{idx2 + 1}-{idx3 + 1}-{idx4 + 1}",
              atom_idx1=idx1,
              atom_idx2=idx2,
              atom_idx3=idx3
            )
  if mol != None:
    res.label = res.get_info(atoms_obj=mol)

  return res



