from dataclasses import dataclass, field
import numpy as np
from typing import Any, Dict, List, Tuple, Union

from ase import Atoms

import ase_utils as au
import features as ft
from features_constraints import Constraint, FeatureConstraint, MoleculeFeatureConstraints
from target_molecule import TargetMoleculeBonds as tm


def validate_structure(
                        mol: Atoms,
                        constraints: List[FeatureConstraint]
                      ) -> bool:

  return MoleculeFeatureConstraints(
              label="",
              molecule=mol,
              feature_constraints=constraints
            ).satisfied()



def target_mol_constraints() -> List[FeatureConstraint]:
  return [
    FeatureConstraint(
        label="",
        feature=ft.Angle("", 0, 0, 0),
        constraint=Constraint("", "<", 0)
      )


  ]

def get_missing_bonds_info(
                            mol: Atoms,
                            bonds: List[Tuple[int, int]]
                          ) -> List[str]:

  res = [f"{x}, dist_in_target={round(tm.bonds_distances[x], 3)}, dist_in_mol={round(mol.get_distance(x[0], x[1]),3)}" for x in bonds]

  return res


def validate_with_target_molecule_bonds(mol: Atoms) -> Dict:
  target_cc = au.convert_bonds_to_set(tm.CC)
  mol_cc = au.get_bonds_as_set(mol=mol, element1="C", element2="C")
  missing_cc_in_target, missing_cc_in_mol = au.get_bonds_diff(bonds1=target_cc, bonds2=mol_cc)
  num_missing_cc = len(missing_cc_in_target) + len(missing_cc_in_mol)

  target_ch = au.convert_bonds_to_set(tm.CH)
  mol_ch = au.get_bonds_as_set(mol=mol, element1="C", element2="H")
  missing_ch_in_target, missing_ch_in_mol = au.get_bonds_diff(bonds1=target_ch, bonds2=mol_ch)
  num_missing_ch = len(missing_ch_in_target) + len(missing_ch_in_mol)

  target_co = au.convert_bonds_to_set(tm.CO)
  mol_co = au.get_bonds_as_set(mol=mol, element1="C", element2="O")
  missing_co_in_target, missing_co_in_mol = au.get_bonds_diff(bonds1=target_co, bonds2=mol_co)
  num_missing_co = len(missing_co_in_target) + len(missing_co_in_mol)

  total_missing_in_target = len(missing_cc_in_target) + len(missing_ch_in_target) + len(missing_co_in_target)
  total_missing_in_mol = len(missing_cc_in_mol) + len(missing_ch_in_mol) + len(missing_co_in_mol)
  is_ok = total_missing_in_target + total_missing_in_mol == 0

  res = { "is_valid": is_ok }

  if num_missing_cc > 0:
    res["C-C"] = {}

  if num_missing_ch > 0:
    res["C-H"] = {}

  if num_missing_co > 0:
    res["C-O"] = {}

  if len(missing_cc_in_target) > 0:
    res["C-C"]["new_cc_bonds"] = get_missing_bonds_info(mol, missing_cc_in_target)

  if len(missing_cc_in_mol) > 0:
    res["C-C"]["lost_cc_bonds"] = get_missing_bonds_info(mol, missing_cc_in_mol)

  if len(missing_ch_in_target) > 0:
    res["C-H"]["new_ch_bonds"] = get_missing_bonds_info(mol, missing_ch_in_target)

  if len(missing_ch_in_mol) > 0:
    res["C-H"]["lost_ch_bonds"] = get_missing_bonds_info(mol, missing_ch_in_mol)

  if len(missing_co_in_target) > 0:
    res["C-O"]["new_co_bonds"] = get_missing_bonds_info(mol, missing_co_in_target)

  if len(missing_co_in_mol) > 0:
    res["C-O"]["lost_co_bonds"] = get_missing_bonds_info(mol, missing_co_in_mol)

  return res



