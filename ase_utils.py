from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
import numpy as np

from ase import Atoms

import constants as C
from dataset import Dataset
import utils as ut
import xyz_parser


def create_ase_atoms_obj(
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


def create_ase_atoms_obj_from_xyz_data(xyz_data: Dict) -> Atoms:

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


  return create_ase_atoms_obj(
                              atomic_nrs=atomic_nrs,
                              coords=elem_coords["coords"],
                              info=info
                            )


def create_ase_atom_objects_from_xyz_file(input_path: Path, name: str) -> List[Atoms]:
  add_info = {
    "name": name,
    }

  res = [create_ase_atoms_obj_from_xyz_data({**x, **add_info}) for x in xyz_parser.read_xyz_file(
                                                              input_path=input_path,
                                                              convert_coords_to_float=True
                                                            )
        ]
  return res


def create_ase_atom_objects_from_dataset(dataset: Dataset) -> List[Atoms]:
  res = []
  for i in range(len(dataset)):
    mols_from_file = create_ase_atom_objects_from_xyz_file(
                                                            input_path=dataset[i],
                                                            name=dataset.names[i]
                                                          )
    res.append(mols_from_file)

  return ut.flatten_list(res)
