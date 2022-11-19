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


