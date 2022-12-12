# Creates dataset from optimisation results

from copy import deepcopy
from datetime import datetime
from itertools import chain
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import ase_utils as au
import constants as C
import utils as ut


def find_row_idx_in_dataframe(
                              df: pd.DataFrame,
                              column: str,
                              value: Any
                            ) -> int:
  try:
    return df[column].loc[lambda x: x==value].index[0]
  except:
    return -1


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


def convert_file_name_to_display_name(file_name: str) -> str:
    res = file_name.lower().replace("mol24_tpsstpss_gd3", "original")
    res = res.replace("mol24_b3p86", "original")
    res = res.replace("mol24_wb97xd", "original")
    res = res.replace("mol24_", "")
    res = res.replace("crest_ex", "ex")
    res = res.replace("crest_", "c")
    res = res.replace("gfn2_", "").replace("gfnff_", "")
    res = res.replace("bp86_", "").replace("pbe1pbe_", "").replace("tpsstpss_gd3_", "")
    res = res.replace("b3p86_", "").replace("wb97xd_", "")
    res = res.replace("def2svpp_svpfit_", "").replace("sto3g_", "").replace("cc_pvtz_tzvpfit_", "")
    res = res.replace("_def2svpp_svpfit", "").replace("_sto3g", "").replace("_cc_pvtz_tzvpfit", "")
    res = res.replace("crest", "c")
    return res


def prepare_dataframe(
                      df: pd.DataFrame,
                      drop_duplicates: bool=True,
                      best_item_idx: int= -1,
                    ) -> pd.DataFrame:

  best_final_energy = df["energy_end"][best_item_idx] \
                        if best_item_idx >= 0 and best_item_idx < len(df) \
                        else min(df["energy_end"])

  res = df.copy()
  res["energy_diff_to_best_au"] = res["energy_end"] - best_final_energy
  res["energy_diff_to_best_kj_mol"] = res["energy_diff_to_best_au"] * C.hartree_in_kJ_per_mol
  res["energy_diff_to_best_kcal_mol"] = res["energy_diff_to_best_au"] * C.hartree_in_kcal_per_mol
  res["source"] = res["input_path"].apply(lambda x: Path(x).name)
  res["label"] = res["input_path"].apply(lambda x: convert_file_name_to_display_name(Path(x).stem))

  if drop_duplicates:
    res.drop_duplicates("label", keep="first", inplace=True)

  res.sort_values("energy_end", inplace=True)
  res.reset_index(inplace=True)

  return res



def find_common_items(
                      list1: List[Any],
                      list2: List[Any]
                    ) -> List[Any]:

  return list(set(list1).intersection(list2))


def prepare_dataframes_by_dft(
                                df: pd.DataFrame,
                                dft_names: List[str]
                              ) -> Dict:

  res = {
      "common_labels": {},
    }
  key_col = "label"
  
  labels = []
  for dft in dft_names:
    mask = df["dft_functional"] == dft
    res[dft] = prepare_dataframe(df=df[mask])
    labels.append(res[dft][key_col])

  for i, dft1 in enumerate(dft_names[:-1]):
    for j in range(i + 1, len(dft_names)):
      dft2 = dft_names[j]
      common_labels = find_common_items(
                        list1=list(res[dft1][key_col]),
                        list2=list(res[dft2][key_col])
                        )
      res["common_labels"][(dft1, dft2)] = common_labels
      res["common_labels"][(dft2, dft1)] = common_labels
          
  return res


def find_common_items_in_dataframes(
                                    column: str,
                                    df_list: List[pd.DataFrame]
                                  ) -> List[Any]:

  '''
    Return the list of common items among all dataframes in df_list.
    Search in the column specified.
    Return unique common items from the column spcified.
  '''

  res = list(set(list(df_list[0][column])))
  for df in df_list[1:]:
    res = find_common_items(
        list1=res,
        list2=list(df[column])
  )

  return res









