from dataclasses import dataclass, field
import numpy as np
from typing import Any, Dict, List, Tuple, Union
from tqdm import tqdm

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
                            bonds: List[Tuple[int, int]],
                            target_distances: Dict,
                          ) -> List[str]:

  symbols = mol.get_chemical_symbols()
  res = []
  for x in bonds:
    item_info = {
      "bond": f"{symbols[x[0]]}-{symbols[x[1]]}",
      "atom_idxs": x,
      "dist_in_target": round(target_distances[x], 3) if x in target_distances else 'n/a',
      "dist_in_mol": round(mol.get_distance(x[0], x[1]),3),
    }
    res.append(item_info)

  return res


def validate_bonds_of_one_mol_with_target_molecule(mol: Atoms) -> Dict:

  target_cc = au.convert_bonds_to_set(tm.CC)
  mol_cc = au.get_bonds_as_set(mol=mol, element1="C", element2="C")
  missing_cc_in_target, missing_cc_in_mol = au.get_bonds_diff(bonds1=target_cc, bonds2=mol_cc)

  target_ch = au.convert_bonds_to_set(tm.CH)
  mol_ch = au.get_bonds_as_set(mol=mol, element1="C", element2="H")
  missing_ch_in_target, missing_ch_in_mol = au.get_bonds_diff(bonds1=target_ch, bonds2=mol_ch)

  target_co = au.convert_bonds_to_set(tm.CO)
  mol_co = au.get_bonds_as_set(mol=mol, element1="C", element2="O")
  missing_co_in_target, missing_co_in_mol = au.get_bonds_diff(bonds1=target_co, bonds2=mol_co)

  missing_in_target = list(missing_cc_in_target) + list(missing_ch_in_target) + list(missing_co_in_target)
  missing_in_mol = list(missing_cc_in_mol) + list(missing_ch_in_mol) + list(missing_co_in_mol)

  num_missing_in_target = len(missing_in_target)
  num_missing_in_mol = len(missing_in_mol)
  is_ok = num_missing_in_target + num_missing_in_mol == 0

  res = {
    "is_valid": is_ok,
    "units": "angstrom",
    "new": [],
    "lost": [],
    }

  if num_missing_in_target > 0:
    res["new"] = get_missing_bonds_info(
                                          mol=mol,
                                          bonds=missing_in_target,
                                          target_distances=tm.bonds_distances
                                        )

  if num_missing_in_mol > 0:
    res["lost"] = get_missing_bonds_info(
                                          mol=mol,
                                          bonds=missing_in_mol,
                                          target_distances=tm.bonds_distances
                                        )


  return res


def validate_bonds_of_many_mols_with_target_molecule(
                  mols_to_validate: List[Atoms],
                  target_mol: Atoms=None
                ) -> List[Dict]:

  progress_bar = tqdm(mols_to_validate)
  invalid_items = []
  #for experiment_dir in tqdm(experimentr_dirs, desc="Processing gaussian log files"):
  for mol in progress_bar:
    progress_bar.set_description(au.get_name_from_atoms(mol))
    validation_result = validate_bonds_of_one_mol_with_target_molecule(mol)
    if not validation_result["is_valid"]:
      validation_result["name"] = au.get_name_from_atoms(mol)
      validation_result["description"] = au.get_description_from_atoms(mol)
      invalid_items.append(validation_result)
  
  summary = {
    "num_mols_total": len(mols_to_validate),
    "num_mols_valid": len(mols_to_validate) - len(invalid_items),
    "num_mols_invalid": len(invalid_items),

  }

  res = {
    "summary": summary,
    "invalid_mols": invalid_items,
  }

  return res

