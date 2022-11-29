import copy
from datetime import datetime
from itertools import chain
import json
import numpy as np
import os
from pathlib import Path
import shutil
from typing import Any, Dict, List, Tuple, Union

################################################################
# Generic functions
################################################################

def copy_files(
                source_folder: Union[str, Path],
                target_folder: Union[str, Path],
                file_extension: str
              ) -> None:

  files_to_copy = get_file_paths_in_dir(
                                        search_dir=Path(source_folder),
                                        file_extension=file_extension
                                        )

  for f in files_to_copy:
    shutil.copyfile(f, Path(target_folder).joinpath(f.name))


def convert_list_of_tuples_to_tuple_of_lists(list_of_tuples: List[Tuple]) -> Tuple[List]:
  return list(zip(*list_of_tuples))


def flatten_list(list_of_lists: List[List[Any]]) -> List[Any]:
  return list(chain.from_iterable(list_of_lists))


def get_file_paths_in_dir(
                          search_dir: Union[str, Path],
                          file_extension: str
                        ) -> List[Path]:

  return [x for x in Path(search_dir).iterdir() \
                if x.is_file() and x.suffix.endswith(file_extension)]


def get_file_paths_in_many_dirs(
                          search_dirs: Union[List[str], List[Path]],
                          file_extension: str
                        ) -> List[Path]:

  return flatten_list(
      [get_file_paths_in_dir(x, file_extension) for x in search_dirs]
    )


def is_within_tolerance(value1: float, value2: float, tolerance: float) -> bool:
  return abs(value1 - value2) < tolerance


def to_json_str(data: Any) -> str:
  return json.dumps(data, sort_keys=False, indent=4)


def write_text_file(file_name: Union[str, Path], text: str) -> str:
  input_path = Path(file_name)
  _ = input_path.write_text(text)
  return str(input_path)


def write_text_file_json(file_name: Union[str, Path], data: Any) -> str:
  return write_text_file(file_name, to_json_str(data))


def read_text_file_as_lines(file_path: Union[str, Path]) -> List[str]:
  with open(file_path) as f:
    lines = f.readlines()

  return lines


def write_text_file_from_lines(
                                file_path: Union[str, Path],
                                lines: List[str]
                                ) -> str:

  with open(file_path, 'w') as f:
    for line in lines:
      if isinstance(line, str):
        f.write(line)
        if len(line) < 2 or ("\n" not in line[-1:]):
          f.write("\n")
      else:
        f.write(f"{line}")
        f.write("\n")

  return str(file_path)


def get_block_start_line_nrs(
                              lines: List[str],
                              search_text: str,
                              ) -> List[int]:

  return [i for i,x in enumerate(lines) if search_text in x]


def get_line_nrs_starts_with_text(
                              lines: List[str],
                              search_text: str,
                              ) -> List[int]:

  return [i for i,x in enumerate(lines) if x.startswith(search_text)]
