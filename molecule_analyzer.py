import copy
from dataclasses import dataclass, field
from datetime import datetime
import json
import numpy as np
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import ase
import ase.io as ase_io
from ase import Atoms
from ase import constraints as ase_constraints
from ase.constraints import FixInternals

import features as ft

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
@dataclass
class SeriesItem():
  name: str
  labels: List[str] = field(default_factory=list)
  values: List[float] = field(default_factory=list)


def calculate_values_groupby_features(
                                      molecules: List[Atoms],
                                      features_list: List[ft.Feature]
                                      ) -> List[SeriesItem]:

  res = []
 
  labels = [mol.info["name"] for mol in molecules]

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

