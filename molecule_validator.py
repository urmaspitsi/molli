from dataclasses import dataclass, field
import numpy as np
from typing import Any, Dict, List, Tuple, Union

from ase import Atoms

import ase_utils as au
import features as ft
from features_constraints import Constraint, FeatureConstraint, MoleculeFeatureConstraints
from target_molecule import TargetMoleculeBonds


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


def validate_with_target_molecule_bonds(mol: Atoms) -> Dict:
  target_cc = au.convert_bonds_to_set(TargetMoleculeBonds.CC)
  mol_cc = au.get_bonds_as_set(mol=mol, element1="C", element2="C")
  missing_cc_in_target, missing_cc_in_mol = au.get_bonds_diff(bonds1=target_cc, bonds2=mol_cc)

  target_ch = au.convert_bonds_to_set(TargetMoleculeBonds.CH)
  mol_ch = au.get_bonds_as_set(mol=mol, element1="C", element2="H")
  missing_ch_in_target, missing_ch_in_mol = au.get_bonds_diff(bonds1=target_ch, bonds2=mol_ch)

  target_co = au.convert_bonds_to_set(TargetMoleculeBonds.CO)
  mol_co = au.get_bonds_as_set(mol=mol, element1="C", element2="O")
  missing_co_in_target, missing_co_in_mol = au.get_bonds_diff(bonds1=target_co, bonds2=mol_co)

  total_missing_in_target = len(missing_cc_in_target) + len(missing_ch_in_target) + len(missing_co_in_target)
  total_missing_in_mol = len(missing_cc_in_mol) + len(missing_ch_in_mol) + len(missing_co_in_mol)
  is_ok = total_missing_in_target + total_missing_in_mol == 0

  res = {
      "is_valid": is_ok,
      "cc": {
              "new_cc_bonds": list(missing_cc_in_target),
              "lost_cc_bonds": list(missing_cc_in_mol),
              },

      "ch": {
              "new_ch_bonds": list(missing_ch_in_target),
              "lost_ch_bonds": list(missing_ch_in_mol),
              },

      "co": {
              "new_co_bonds": list(missing_co_in_target),
              "lost_co_bonds": list(missing_co_in_mol),
              },
    }

  return res



