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


def write_xyz_file_from_dict(
                              output_path: Union[str, Path],
                              xyz_as_list_of_dicts: List[Dict]
                            ) -> int:
  '''
    Writes xyz dictionary into text file.
    Reverse function to 'read_xyz_file(...)'
  '''

  pass


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
            }
        )

      start_line = end_line
      idx +=1

    else: # empty or non-numeric line
      start_line += 1

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

