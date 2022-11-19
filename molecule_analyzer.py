import copy
from dataclasses import dataclass
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
#

