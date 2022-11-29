from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
import numpy as np

from ase import Atoms
from ase.build import minimize_rotation_and_translation

import constants as C
from dataset import Dataset
import utils as ut
import xyz_parser


def align_2_molecules_min_rmsd(target: Atoms, atoms_to_align: Atoms) -> Atoms:
  '''
    Aligns atoms_to_align with target, based on min rmsd.
    Returns copy of atoms_to_align with new position coordinates.
  '''

  res = deepcopy(atoms_to_align)
  minimize_rotation_and_translation(target, res)

  return res


def create_ase_atoms(
                      atomic_nrs: List[int],
                      coords: List[List[float]],
                      info: Dict,
                    ) -> Atoms:

  # TODO: initial atoms with cell definition and pbc: https://wiki.fysik.dtu.dk/ase/ase/atoms.html
  atoms_dict = {
      "numbers": atomic_nrs,
      "positions": coords,
      "cell": np.zeros((3, 3), dtype=float), # .fill(1.0)
      "pbc": np.zeros((3,), dtype=float),    # .fill(1)
    }

  res = Atoms.fromdict(atoms_dict)

  if info != None and len(info) > 0:
    if res.info != None and len(res.info) > 0:
      res.info = { **res.info, **info }
    else:
      res.info = info

  return res


def create_ase_atoms_from_xyz_data(xyz_data: Dict) -> Atoms:

  xyz_lines_are_strings = True \
                            if isinstance(xyz_data["xyz_lines"][0], str) \
                            else False

  elem_coords = xyz_parser.split_xyz_lines_elements_coords_from_str(xyz_data["xyz_lines"]) \
                  if xyz_lines_are_strings \
                  else xyz_parser.split_xyz_lines_elements_coords_from_tuple(xyz_data["xyz_lines"])

  atomic_nrs = [C.get_atomic_number(i) for i in elem_coords["elements"]]
  info = {}
  info_keys = ["name", "description"]

  for k in info_keys:
    if k in xyz_data:
      info[k] = xyz_data[k]

  return create_ase_atoms(
                          atomic_nrs=atomic_nrs,
                          coords=elem_coords["coords"],
                          info=info
                        )


def create_ase_atoms_list_from_xyz_file(
                                          input_path: Path,
                                          name: str
                                        ) -> List[Atoms]:

  parsed_xyz_data = [x for x in xyz_parser.read_xyz_file(
                                              input_path=input_path,
                                              convert_coords_to_float=True
                                            )
        ]
  add_item_nr_to_name = True if len(parsed_xyz_data) > 1 else False

  # multiple xyz blocks in the file: trajectory, optimization steps, conformers etc.
  # append idx+1 numeric value to the name: idx=0, "name" -> "name_1"
  # idx is zero based, added number will be 1-based.
  for xyz in parsed_xyz_data:
    if add_item_nr_to_name and "idx" in xyz:
      item_nr = xyz["idx"] + 1
      xyz["name"] = f"{name}_{item_nr}"
    else:
      # single xyz geometry or no idx information
      xyz["name"] = name

  res = [create_ase_atoms_from_xyz_data(x) for x in parsed_xyz_data]

  return res


def create_ase_atoms_list_from_dataset(dataset: Dataset) -> List[Atoms]:

  res = []
  for i in range(len(dataset)):
    mols_from_file = create_ase_atoms_list_from_xyz_file(
                                                          input_path=dataset[i],
                                                          name=dataset.names[i]
                                                        )
    res.append(mols_from_file)

  return ut.flatten_list(res)


def compare_if_molecules_are_equal(mol1: Atoms, mol2: Atoms) -> bool:
  return len(mol1.numbers) == len(mol2.numbers) \
        and (mol1.numbers == mol2.numbers).all() \
        and np.allclose(
              mol1.get_all_distances(),
              mol2.get_all_distances(),
              rtol=1e-05,
              atol=1e-08,
              equal_nan=True
            )


def rmsd_of_distances(mol1: Atoms, mol2: Atoms) -> float:
  d1 = mol1.get_all_distances()
  d2 = mol2.get_all_distances()
  return np.sqrt(np.mean((d1 - d2)**2))


def try_get_info_item_from_atoms(mol: Atoms, info_key: str) -> str:
  try:
    if info_key in mol.info:
      return mol.info[info_key]
    else:
      return ""
  except:
    return ""


def try_get_name_from_atoms(mol: Atoms) -> str:
  return try_get_info_item_from_atoms(mol, info_key="name")


def try_get_description_from_atoms(mol: Atoms) -> str:
  return try_get_info_item_from_atoms(mol, info_key="description")


def write_ase_atoms_to_xyz_file(
            atoms_list: List[Atoms],
            output_path: Union[str, Path]
            ) -> None:

  res = []
  for mol in atoms_list:
    description = f"{try_get_name_from_atoms(mol)} {try_get_description_from_atoms(mol)}"
    lines = [xyz_parser.convert_xyz_coords_to_str(el, x, y, z) \
              for el, (x,y,z) in zip(mol.get_chemical_symbols(), list(mol.positions))]

    res.append(len(lines))
    res.append(description)
    res.extend(lines)

  return ut.write_text_file_from_lines(file_path=output_path, lines=res)
