from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

from dataset import Dataset


mol24_crest_conformers = Dataset.create_with_names(
  description="mol24_crest_conformers",
  paths_with_names=
      [
        (Path("C:/tmp/Molli/data/mol24_final.xyz"), "mol24_final"),
        (Path("C:/tmp/Molli/data/mol24_ex5_gfnff_crest_conformers.xyz"), "mol24_ex5_gfnff_crest_conformers"),
        (Path("C:/tmp/Molli/data/mol24_ex8_gfn2_crest_conformers.xyz"), "mol24_ex8_gfn2_crest_conformers"),
        (Path("C:/tmp/Molli/data/mol24_ex10_gfnff_crest_conformers.xyz"), "mol24_ex10_gfnff_crest_conformers"),
        (Path("C:/tmp/Molli/data/mol24_ex15_gfnff_crest_conformers.xyz"), "mol24_ex15_gfnff_crest_conformers"),
        (Path("C:/tmp/Molli/data/mol24_ex16_gfn2_crest_conformers.xyz"), "mol24_ex16_gfn2_crest_conformers"),
      ]
)


mol24_xtb_gfn2_opt_steps = Dataset.create_with_names(
  description="mol24_xtb_gfn2_optimize",
  paths_with_names=
      [
        (Path("C:/tmp/Molli/data/mol24_ex17_gfn2_opt_steps.xyz"), "mol24_ex17_gfn2_opt_steps"),
      ]
)


mol24_xtb_gfnff_opt_steps = Dataset.create_with_names(
  description="mol24_xtb_gfnff_optimize",
  paths_with_names=
      [
        (Path("C:/tmp/Molli/data/mol24_ex18_gfnff_opt_steps.xyz"), "mol24_ex18_gfnff_opt_steps"),
      ]
)

