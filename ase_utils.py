from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, Dict, List, Set, Tuple, Union
import numpy as np

from ase import Atoms
from ase.build import minimize_rotation_and_translation
from ase.geometry.analysis import Analysis

import constants as C
from dataset import Dataset
from file_source import MultiItemFileSource
import utils as ut
import metrics as ms
import xyz_parser


def align_2_molecules_min_rmsd(
                                target: Atoms,
                                atoms_to_align: Atoms,
                              ) -> Atoms:

  '''
    Aligns atoms_to_align with target, based on min rmsd.
    Returns copy of atoms_to_align with new position coordinates.
  '''

  res = deepcopy(atoms_to_align)
  minimize_rotation_and_translation(target, res)

  return res


def calculate_metric_between_two_molecules(
                    target: Atoms,
                    mol: Atoms,
                    align: bool,
                    metric_function: Callable,
                  ) -> float:
  '''
    Applies metric_function on target and mol:
      metric_function(mol1=target, mol2=mol) -> float
    metric_function must:
      take 2 input parameters: mol1: Atoms, mol2: Atoms
      and return float.
    
    align: True/False, whether to align target and mol before metric calculation.
    alignment is done by align_2_molecules_min_rmsd()
  '''

  aligned_mol = align_2_molecules_min_rmsd(
                                            target=target,
                                            atoms_to_align=mol
                                          ) if align else mol

  return metric_function(
                          mol1=target,
                          mol2=aligned_mol
                        )


def calculate_metric_one_to_many(
                                target: Atoms,
                                mols: List[Atoms],
                                align: bool,
                                metric_function: Callable,
                              ) -> List[float]:
  '''
    Applies metric_function on target and mol:
      metric_function(mol1=target, mol2=mol) -> float
    metric_function must:
      take 2 input parameters: mol1: Atoms, mol2: Atoms
      and return float.
    
    align: True/False, whether to align target and mol before metric calculation.
    alignment is done by align_2_molecules_min_rmsd()
  '''

  return [calculate_metric_between_two_molecules(
                target=target,
                mol=m,
                align=align,
                metric_function=metric_function
                ) for m in mols]


def calculate_metric_many_to_many(
                                  targets: List[Atoms],
                                  mols: List[Atoms],
                                  align: bool,
                                  metric_function: Callable,
                                ) -> List[List[float]]:
  '''
    Applies metric_function on target and mol:
      metric_function(mol1=target, mol2=mol) -> float
    metric_function must:
      take 2 input parameters: mol1: Atoms, mol2: Atoms
      and return float.
    
    align: True/False, whether to align target and mol before metric calculation.
    alignment is done by align_2_molecules_min_rmsd()
  '''

  return [calculate_metric_one_to_many(
            target=t,
            mols=mols,
            align=align,
            metric_function=metric_function
            ) for t in targets]


def calculate_metric_cross(
                            mols: List[Atoms],
                            align: bool,
                            metric_function: Callable,
                          ) -> Dict:
  '''
    Applies metric_function on every item in mols with every other item in mols:
      metric_function(mol1=mols[i], mol2=mols[j]) -> float
      i,j = 0 to len(mols)
    metric_function must:
      take 2 input parameters: mol1: Atoms, mol2: Atoms
      and return float.
    
    align: True/False, whether to align target and mol before metric calculation.
    alignment is done by align_2_molecules_min_rmsd()

    returns Dictionary {
      names: List[str],
      values: Dict[Tuple[int,int], float]
    }
       
  '''
  res = {
    "info": [
              {
                "name": get_name_from_atoms(mol=mol),
                "description": get_description_from_atoms(mol=mol),
                "source": get_source_from_atoms(mol=mol),
              } for mol in mols
            ],

    "values": {},
  }

  for i, mol1 in enumerate(mols[:-1]):
    for j in range(i + 1, len(mols)):
      res["values"][(i ,j)] = calculate_metric_between_two_molecules(
        target=mol1,
        mol=mols[j],
        align=align,
        metric_function=metric_function
      )

  return res


def calculate_metric_xyz_files(
                              target_mols_path: Path,
                              mols_path: Path,
                              align: bool,
                              metric_function: Callable,
                            ) -> List[List[float]]:

  target_mols = create_ase_atoms_list_from_xyz_file(
                                                    input_path=target_mols_path,
                                                    name=target_mols_path.stem
                                                  )

  mols = create_ase_atoms_list_from_xyz_file(
                                              input_path=mols_path,
                                              name=mols_path.stem
                                            )

  return calculate_metric_many_to_many(
                                    targets=target_mols,
                                    mols=mols,
                                    align=align,
                                    metric_function=metric_function
                                  )


def calculate_metric_between_xyz_files(
                                      target_xyz_path: Path,
                                      xyz_path: Path,
                                      max_value: float,
                                      align: bool,
                                      metric_function: Callable,
                                    ) -> Dict:

  values = calculate_metric_xyz_files(
                      target_mols_path=target_xyz_path,
                      mols_path=xyz_path,
                      align=align,
                      metric_function=metric_function
                    )

  less_than_max_value = []
  is_same_file = True if target_xyz_path == xyz_path else False
  for i, vals_i in enumerate(values):
    j_start = i + 1 if is_same_file else 0
    for j in range(j_start, len(vals_i)):
      rmsd_val = values[i][j]
      if rmsd_val < max_value:
        less_than_max_value.append(("target_source_value", i, j, rmsd_val))

  less_than_max_value = sorted(less_than_max_value, key=lambda x: x[3])

  res = {
    "less_than_max_value": less_than_max_value,
    "all_values": values,
  }

  return res


def calculate_metric_xyz_file(
                              xyz_path: Path,
                              max_value: float,
                              align: bool,
                              metric_function: Callable
                            ) -> Dict:

  res = calculate_metric_between_xyz_files(
      target_xyz_path=xyz_path,
      xyz_path=xyz_path,
      max_value=max_value,
      align=align,
      metric_function=metric_function
    )

  keys = set()
  less_than_max_value = []
  for name,i,j,val in res["less_than_max_value"]:
    if not ((i, j) in keys or (j, i)) in keys:
      less_than_max_value.append((name, i, j, val))
      keys.add((i, j))

  res["less_than_max_value"] = less_than_max_value

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
  info_keys = ["name", "description", "source"]

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
                                          name: str,
                                          item_idxs: Union[List[int], int, None]=None
                                        ) -> List[Atoms]:
  '''
    Returns list of Atoms from the xyz file specified.
    If item_idxs is of type List[int], then only these items are returned: 0-based list indices.
    If item_idxs is of type int, then returns number of items eaqual to the value of item_idxs.
    First item is 0th element and last item is the last item in the file.
    All intermediate items are linearly spaced.
    0-index item, ..., ..., last-item in the file
    Idea is to get a sub-sequence of certain lenght from a list,
    starting from the first item and ending with the last item in the original list.
    Useful for generating equal length lists when comparing different opt trajectories.
  '''

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

  idxs = None
  # If item_idxs is integer then means num_items and generates linspaced 
  # list of integers between 0 and len(parsed_xyz_data)
  if item_idxs:
    if isinstance(item_idxs, int) and item_idxs > 1:
      idxs = ut.linspace_idxs(start_idx=0, end_idx=(len(parsed_xyz_data) - 1), num_items=item_idxs)
    elif isinstance(item_idxs, list) and len(item_idxs) > 0:
      idxs = item_idxs

  # Filters by given list of item_idxs
  if idxs and len(idxs) > 0:
    parsed_xyz_data = ut.get_list_slice_by_idxs(
                          input_list=parsed_xyz_data,
                          idxs=idxs
                          )

  res = [create_ase_atoms_from_xyz_data(x) for x in parsed_xyz_data]

  return res


def create_ase_atoms_list_from_xyz_files(
                                          input_paths: List[Path]
                                        ) -> List[Atoms]:

  return ut.flatten_list(
              [
                create_ase_atoms_list_from_xyz_file(
                                                    input_path=p,
                                                    name=p.stem
                                                  )
                    for p in input_paths
              ]
            )


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


def extract_mols_from_xyz_files(
      sources: List[MultiItemFileSource]
    ) -> List[List[Atoms]]:

  res = [create_ase_atoms_list_from_xyz_file(
                    input_path=src.file_path,
                    name=src.name,
                    item_idxs=src.item_idxs)
          for src in sources]

  return res


def get_info_item_from_atoms(mol: Atoms, info_key: str) -> str:
  try:
    if info_key in mol.info:
      return mol.info[info_key]
    else:
      return ""
  except:
    return ""


def get_name_from_atoms(mol: Atoms) -> str:
  return get_info_item_from_atoms(mol, info_key="name")


def get_description_from_atoms(mol: Atoms) -> str:
  return get_info_item_from_atoms(mol, info_key="description")


def get_source_from_atoms(mol: Atoms) -> str:
  return get_info_item_from_atoms(mol, info_key="source")


def write_ase_atoms_to_xyz_file(
            atoms_list: List[Atoms],
            output_path: Union[str, Path]
            ) -> str:

  res = []
  for mol in atoms_list:
    description = f"{get_name_from_atoms(mol)}, {get_description_from_atoms(mol)}, source: {get_source_from_atoms(mol)}"
    lines = [xyz_parser.convert_xyz_coords_to_str(el, x, y, z) \
              for el, (x,y,z) in zip(mol.get_chemical_symbols(), list(mol.positions))]

    res.append(len(lines))
    res.append(description)
    res.extend(lines)

  return ut.write_text_file_from_lines(file_path=output_path, lines=res)


def write_aligned_xyz_file(
            input_path: Union[str, Path],
            output_path: Union[str, Path],
          ):

  mols = create_ase_atoms_list_from_xyz_file(
                  input_path=input_path,
                  name=input_path.stem
                )

  aligned_mols = [
      align_2_molecules_min_rmsd(
                                  target=mols[0],
                                  atoms_to_align=x
                                ) for x in mols
    ]

  return write_ase_atoms_to_xyz_file(
            atoms_list=aligned_mols,
            output_path=output_path
          )


def get_bonds(
              mol: Atoms,
              element1: str,
              element2: str,
            ) -> List[Tuple[int]]:
  '''
    Returns list of tuples of atom idxs.
  '''

  ana = Analysis(mol)
  return ana.get_bonds(element1, element2, unique=True)[0]


def get_bonds_as_set(
              mol: Atoms,
              element1: str,
              element2: str,
            ) -> List[Tuple[int]]:
  '''
    Returns list of tuples of atom idxs.
  '''

  return convert_bonds_to_set(
                              bonds=get_bonds(
                                      mol=mol,
                                      element1=element1,
                                      element2=element2
                                    )
                            )


def convert_bonds_to_set(bonds: List[Tuple[int]]) -> Set[Tuple[int]]:
  res = set()
  for i, j in bonds:
    if i < j:
      res.add((i, j))
    else:
      res.add((j, i))

  return res


def get_bonds_diff(bonds1: Set[Tuple[int]], bonds2: Set[Tuple[int]]) -> Tuple[Set[Tuple[int]]]:
  '''
    Returns tuple of 2-sets: missing_in_bonds1, missing_in_bonds2
  '''
  diff = bonds1.symmetric_difference(bonds2)
  missing_in_bonds1 = set()
  missing_in_bonds2 = set()
  for pair in diff:
    if pair not in bonds1:
      missing_in_bonds1.add(pair)
    else:
      missing_in_bonds2.add(pair)

  return (missing_in_bonds1, missing_in_bonds2)


