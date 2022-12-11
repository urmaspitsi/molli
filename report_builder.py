# Creates dataset from optimisation results

from copy import deepcopy
from datetime import datetime
from itertools import chain
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import ase_utils as au
import utils as ut


def read_aggregate_summary_to_dict(
          input_path: Union[str, Path]
        ) -> Dict:

  '''
    Reads aggregate summary that is generated from gaussian log files.
    Return Pandas DataFrame object.
  '''    

  data = ut.read_json_file_to_dict(file_path=input_path)

  columns_from_experiment = [
    "input_path",
    ]

  columns_from_results = [
    "final_xyz",
    "opt_steps",
    ]

  columns_from_scf_summary = [
    "gaussian_version",
    "gaussian_command",
    "dft_functional",
    "dft_info",
    "optimization_converged",
    "elapsed_time_str",
    "elapsed_time_minutes",
    "num_steps",
    "minutes_per_step",
    "energy_start",
    "energy_end",
    "energy_delta",
    "energy_delta_text",
    "job_cpu_time",
    "job_cpu_hours",
    "job_cpu_hours_per_step",
    "job_completion_datetime",
    "energy_diff_to_best",
    ]

  columns = columns_from_scf_summary + columns_from_experiment + columns_from_results

  res = {col: [] for col in columns }

  for ex in data["experiments"]:
    for col in columns_from_experiment:
      if col in ex:
        res[col].append(ex[col])
      else:
        res[col].append(None)

    for col in columns_from_results:
      if col in ex["results"]:
        res[col].append(ex["results"][col])
      else:
        res[col].append(None)

    for col in columns_from_scf_summary:
      if col in ex["results"]["scf_summary"]:
        res[col].append(ex["results"]["scf_summary"][col])
      else:
        res[col].append(None)

  return res


def read_aggregate_summary_to_dataframe(
          input_path: Union[str, Path]
        ) -> pd.DataFrame:
  '''
    Reads aggregate summary that is generated from gaussian log files.
    Return Pandas DataFrame object.
  '''    
  return pd.DataFrame(read_aggregate_summary_to_dict(input_path=input_path))



