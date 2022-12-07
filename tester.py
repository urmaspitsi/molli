
# %%

from datetime import datetime
from itertools import chain
import json
import numpy as np
import os
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import ase_utils as au
import metrics as ms
import report_builder as rb
import utils as ut
import xyz_parser

# %%
input_path = Path("C:/tmp/gaussian/aggregated_results/all_summaries/all_summaries_log.txt")
data_dict = rb.read_aggregate_summary_to_dict(input_path=input_path)
df = pd.DataFrame(data_dict)

# %%
df
# %%

# %%
a = np.array([[0,0,0],[1,1,1],[2,2,0]])
b = np.array([[3,4,0],[0,0,0],[3,3,0]])

# %%
b
# %%
np.linalg.norm(a - b, axis=1)
# %%
