import numpy as np
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union
import pandas as pd

from ase import Atoms

import ase_utils as au
from dataset import Dataset
import features as ft
from series_item import SeriesItem

# Dataset: List[Path]
# Add: file path, List[file path], directory path (incl. sub dirs), dir path list
#   generates internally unique List[file path]
# Remove: path

# DataLoader
# Define features/constraints: bonds[a1,a2], angles[a1,a2,a3]
# Filter by features/constraints
# Load single xyz, multi xyz (opt steps xyz)
# Load to memory, filter by items

# Analyzer
# Calculate interesting features/properties
#   output: pandas, csv, dict/json
#   visualize: graph, charts, histograms



# Tasks:
# 1. optimization steps: input 1 xyz file, contains many xyz-blocks (geometries)
#   Q: questions to ask? timeseries of a bond length(s), angle(s).
#    How do these evolve through optimization process?
#   R: report: timeseries plots, overlay plots for angles(?)
#
# 2. Conformers ensembles: input 1 or more xyz files, contain many xyz blocks.
#     each block is independent of the others.
#   Q: questions to ask? filter/group by specific or aggregate features. bonds, angles, etc.
#   R: report: histogram, 
#
#
# Create template datasets: 
#   opt-steps: xtb, gaussian, turbomole: compare various model chemistries
#   conformers: crest, crest-reoptimized
#
#
# Molecules and Features pairs: (mol, feature)?
# Feature: over many molecules: for m in mols: calc_feature_value
#
#
#
#

def analyze_by_features(
                        molecules: List[Atoms],
                        features_list: List[ft.Feature]
                        ) -> Dict:

  res = {
    "info": {
      "num_molecules": len(molecules),
      "num_features": len(features_list),
      "mol_names": [au.get_name_from_atoms(mol=mol) for mol in molecules],
      "features_labels": [ft.label for ft in features_list],

    },
    "results": {}
  }
 
  for ft in features_list:
    res["results"][ft.label] = [
        {
          "name": au.get_name_from_atoms(mol=mol),
          "description": au.get_description_from_atoms(mol=mol),
          "source": au.get_source_from_atoms(mol=mol),
          "value": ft.calculate_value(atoms_obj=mol),
        } for mol in molecules
      ]

  return res


def analyze_by_features_to_dataframe(
                        molecules: List[Atoms],
                        features_list: List[ft.Feature]
                        ) -> pd.DataFrame:

  analysis_res = analyze_by_features(
                    molecules=molecules,
                    features_list=features_list
                  )

  res = pd.DataFrame()
  res["mol_names"] = analysis_res["info"]["mol_names"]

  columns = list(analysis_res["results"].keys())
  for col in columns:
    res[col] = [x["value"] for x in analysis_res["results"][col]]

  return res


def calculate_values_groupby_features(
                                      molecules: List[Atoms],
                                      features_list: List[ft.Feature]
                                      ) -> List[SeriesItem]:

  res = []
 
  labels = [au.get_name_from_atoms(mol=mol) for mol in molecules]

  for ft in features_list:
    features_values = [ft.calculate_value(mol) for mol in molecules]

    series_result = SeriesItem(
                              name=ft.label,
                              labels=labels,
                              values=features_values
                              )

    res.append(series_result)

  return res


def calculate_values_groupby_molecules(
                                      molecules: List[Atoms],
                                      features_list: List[ft.Feature]
                                      ) -> List[SeriesItem]:

  res = []
 
  labels = [ft.label for ft in features_list]

  for mol in molecules:
    features_values = [ft.calculate_value(mol) for ft in features_list]

    series_result = SeriesItem(
                              name=mol.info["name"],
                              labels=labels,
                              values=features_values
                              )

    res.append(series_result)

  return res


def calculate_dataset_list(
                            dataset_list: List[Dataset],
                            features_list: List[ft.Feature]
                        ) -> List[Dict]:

  res = []

  for dset in dataset_list:
    dset_dict = {}
    dset_dict["dataset_description"] = dset.description
    dset_dict["dataset_sources"] = dset.names
    molecules = au.create_ase_atoms_list_from_dataset(dset)
    dset_dict["molecules"] = molecules
    dset_dict["molecules_names"] = [x.info["name"] for x in molecules]

    dset_dict["features_list"] = features_list
    dset_dict["features_info"] = []
    if len(molecules) > 0:
      dset_dict["features_info"] = [x.get_info(molecules[0]) for x in features_list]

    dset_dict["calc_by_features"] = calculate_values_groupby_features(
                                          molecules=molecules,
                                          features_list=features_list
                                        )
    res.append(dset_dict)

  return res
