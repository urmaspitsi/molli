import copy
from datetime import datetime
from itertools import chain
import json
import numpy as np
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
import utils as ut

import constants as C


def convert_xyz_coords_to_str(
                                element: str,
                                x: float,
                                y: float,
                                z: float
                              ) -> str:
  precision = 8
  x_str = str(round(x, precision))
  y_str = str(round(y, precision))
  z_str = str(round(z, precision))

  len_el = len(element)
  len_x = len(x_str)
  len_y = len(y_str)
  len_z = len(z_str)

  space_1 = " " * (20 - len_el - len_x)
  space_2 = " " * (38 - len_el - len(space_1) - len_x - len_y)
  space_3 = " " * (56 - len_el - len(space_1) - len_x - len(space_2) - len_y - len_z)

  return f"{element}{space_1}{x_str}{space_2}{y_str}{space_3}{z_str}"


def convert_xyz_lines_to_dict(
                              xyz_data: List[str],
                              convert_coords_to_float: bool=False,
                            ) -> Dict:

  xyz_lines = [convert_xyz_str_to_coords(x) for x in xyz_data[2:]] \
              if convert_coords_to_float \
              else xyz_data[2:]

  return {
          "num_atoms": xyz_data[0],
          "description": xyz_data[1],
          "xyz_lines": xyz_lines
        }


def convert_xyz_str_to_coords(xyz_row: str) -> Tuple:

  element, x, y, z = xyz_row.strip().split()[:4]
  element = element.strip()
  x = float(x.strip())
  y = float(y.strip())
  z = float(z.strip())

  return (element, [x, y, z])


def read_xyz_file(
                  input_path: Union[str, Path],
                  convert_coords_to_float: bool=False,
                  ) -> List[Dict]:
  ''' Read xyz file into list of dictionaries.
      Returns List[Dict]
      Each dictionary corresponds to 1 xyz block.
      In case xyz file contains 1 xyz block ie file with single geometry,
      then returns list containing 1 element.

      Result dictionary contains keys:
          'num_atoms', 'description', 'xyz_lines', 'start_line'

      convert_coords_to_float=False: xyz-lines remain as strings.
          xyz_lines type: List[str]
      
      convert_coords_to_float=True: xyz-lines are converted into 2-Tuple:
          xyz_lines type: List[Tuple[str, 3-List[float]]]
          Tuple[str, 3-List[float]] is -> (element, [x,y,z]).
  '''

  res = []
  lines = ut.read_text_file_as_lines(file_path=input_path)
  num_lines = len(lines)
  start_line = 0
  idx = 0
  while start_line < num_lines:
    first_line = lines[start_line].strip()
    # first line must be numeric: number of atoms
    if first_line.isdigit():
      num_atoms = int(lines[start_line].strip())
      description = lines[start_line + 1].strip()
      end_line = start_line + num_atoms + 2

      if convert_coords_to_float:
        xyz_lines = [convert_xyz_str_to_coords(x) for x in lines[(start_line + 2):end_line]]
      else:
        xyz_lines = [x.lstrip() for x in lines[(start_line + 2):end_line]]

      res.append(
          {
            "num_atoms": num_atoms,
            "description": description,
            "xyz_lines": xyz_lines,
            "start_line": start_line,
            "idx": idx,
            "source": Path(input_path).name,
            }
        )

      start_line = end_line
      idx +=1

    else: # empty or non-numeric line
      start_line += 1

  return res


def read_xyz_many_files(
                  list_of_input_path_and_idx_tuples: List[Tuple[Union[str, Path], Union[int, List[int]]]],
                  convert_coords_to_float: bool=False,
                  ) -> List[Dict]:
  '''
    Reads input xyz files and aggregates them into one xyz file.
    If idx is specified, then takes only xyz block with that index,
    file contains more than 1 xyz blocks, eg optimization trajectory file,
    or conformers file.

  ''' 
  res = []
  for p, idx in list_of_input_path_and_idx_tuples:
    xyz_file_data = read_xyz_file(
                      input_path=p,
                      convert_coords_to_float=convert_coords_to_float
                      )

    num_xyz_blocks = len(xyz_file_data)

    if isinstance(idx, int) and abs(idx) < num_xyz_blocks:
      res.append(xyz_file_data[idx])
    elif isinstance(idx, list) or isinstance(idx, tuple):
      for i in idx:
        if isinstance(i, int) and abs(i) < num_xyz_blocks:
          res.append(xyz_file_data[i])
    else:
      res.extend(xyz_file_data)

  return res


def split_xyz_lines_elements_coords_from_str(xyz_lines: List[str]) -> Dict:

  element_and_coords = [convert_xyz_str_to_coords(x) for x in xyz_lines]
  elements = [x for x, _ in element_and_coords]
  coords = [x for _, x in element_and_coords]

  return {
          "elements": elements,
          "coords": coords
        }


def split_xyz_lines_elements_coords_from_tuple(xyz_lines: List[Tuple[str, List[float]]]) -> Dict:

  elements, coords = list(zip(*xyz_lines))

  return {
          "elements": elements,
          "coords": coords
        }


def write_xyz_file_from_list_of_dicts(
                              output_path: Union[str, Path],
                              xyz_as_list_of_dicts: List[Dict]
                            ) -> str:
  '''
    Writes xyz dictionary into text file.
    Reverse function to 'read_xyz_file(...)'
  '''
  res = []
  for xyz_dict in xyz_as_list_of_dicts:
    res.append(xyz_dict["num_atoms"])
    res.append(xyz_dict["description"])
    res.extend(xyz_dict["xyz_lines"])

  return ut.write_text_file_from_lines(
          file_path=output_path,
          lines=res
        )


def write_xyz_files_into_one(
      list_of_input_path_and_idx_tuples: List[Tuple[Union[str, Path], Union[int, List[int]]]],
      output_path: Union[str, Path],
    ) -> str:
  '''
    Reads input xyz files and aggregates them into one xyz file.
    If idx is specified, then takes only xyz block with that index,
    file contains more than 1 xyz blocks, eg optimization trajectory file,
    or conformers file.
    Writes the resulting data into one xyz file specified by output_path.
    Returns full path of the output file, if successful.
    Otherwise returns error message as string.
  ''' 

  xyz_data = read_xyz_many_files(
    list_of_input_path_and_idx_tuples=list_of_input_path_and_idx_tuples,
    convert_coords_to_float=False
    )

  return write_xyz_file_from_list_of_dicts(
          output_path=output_path,
          xyz_as_list_of_dicts=xyz_data
          )

