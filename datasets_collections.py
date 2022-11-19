from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

from dataset import Dataset

#-----------------------------------------------------------------------------
# mol24 final: Original input (Öeren 2011)
#-----------------------------------------------------------------------------

mol24_final = Dataset.create_with_names(
  description="mol24_final original input from Öeren 2011",
  paths_with_names=
      [
        (Path("C:/tmp/Molli/data/mol24_final.xyz"), "mol24_final"),
      ]
)


#-----------------------------------------------------------------------------
# mol24 CREST conformers
#-----------------------------------------------------------------------------

mol24_crest_conformers = Dataset.create_with_names(
  description="mol24_crest_conformers",
  paths_with_names=
      [
        (Path("C:/tmp/Molli/data/mol24_ex5_gfnff_crest_conformers.xyz"), "mol24_ex5_gfnff_crest_conformers"),
        (Path("C:/tmp/Molli/data/mol24_ex8_gfn2_crest_conformers.xyz"), "mol24_ex8_gfn2_crest_conformers"),
        (Path("C:/tmp/Molli/data/mol24_ex10_gfnff_crest_conformers.xyz"), "mol24_ex10_gfnff_crest_conformers"),
        (Path("C:/tmp/Molli/data/mol24_ex15_gfnff_crest_conformers.xyz"), "mol24_ex15_gfnff_crest_conformers"),
        (Path("C:/tmp/Molli/data/mol24_ex16_gfn2_crest_conformers.xyz"), "mol24_ex16_gfn2_crest_conformers"),
      ]
)


#-----------------------------------------------------------------------------
# mol24 xtb optimization steps on original mol24_final.xyz
#-----------------------------------------------------------------------------

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


#-----------------------------------------------------------------------------
# CREST best reoptimize: compare first and last step, different model chemistries
#-----------------------------------------------------------------------------

# bp86_def2svpp_svpfit, crest best 1-7
names = [f"mol24_ex16_gfn2_bp86_def2svpp_svpfit_crest_{i}_opt_steps" for i in [1,2,3,4,5,6,7]]

mol24_crest_reopt_gfn2_bp86_def2svpp_svpfit = Dataset.create_with_names(
  description="mol24_crest_reoptimize_ex16_gfn2_bp86_def2svpp_svpfit",
  paths_with_names=[(Path(f"C:/tmp/Molli/data/{n}.xyz"), n) for n in names],
)

# bp86_def2svpp_svpfit, crest best 1-11
names = [f"mol24_ex15_gfnff_bp86_def2svpp_svpfit_crest_{i}_opt_steps" for i in range(1, 12)]

mol24_crest_reopt_gfnff_bp86_def2svpp_svpfit = Dataset.create_with_names(
  description="mol24_crest_reoptimize_ex15_gfnff_bp86_def2svpp_svpfit",
  paths_with_names=[(Path(f"C:/tmp/Molli/data/{n}.xyz"), n) for n in names],
)

# bp86_sto3g, crest best 1-7
names = [f"mol24_ex16_gfn2_bp86_sto3g_crest_{i}_opt_steps" for i in [1,2,3,4,5,6,7]]
          
mol24_crest_reopt_gfn2_bp86_sto3g = Dataset.create_with_names(
  description="mol24_crest_reoptimize_ex16_gfn2_bp86_sto3g",
  paths_with_names=[(Path(f"C:/tmp/Molli/data/{n}.xyz"), n) for n in names],
)

# bp86_sto3g, crest best 1-20
names = [f"mol24_ex15_gfnff_bp86_sto3g_crest_{i}_opt_steps" for i in range(1, 21)]
          
mol24_crest_reopt_gfnff_bp86_sto3g = Dataset.create_with_names(
  description="mol24_crest_reoptimize_ex15_gfnff_bp86_sto3g",
  paths_with_names=[(Path(f"C:/tmp/Molli/data/{n}.xyz"), n) for n in names],
)

# pbe1pbe_sto3g, crest best
mol24_crest_reopt_gfn2_pbe1pbe_sto3g = Dataset.create_with_names(
  description="mol24_crest_reoptimize_ex8_gfn2_pbe1pbe_sto3g",
  paths_with_names=
      [
        (Path("C:/tmp/Molli/data/mol24_ex8_gfn2_crest_best_PBE1PBE_STO3G_opt_opt_steps.xyz"), "mol24_ex8_gfn2_crest_best_PBE1PBE_STO3G_opt_opt_steps"),
      ]
)

# bp86_sto3g, crest best
mol24_crest_reopt_gfn2_bp86_sto3g = Dataset.create_with_names(
  description="mol24_crest_reoptimize_ex8_gfn2_bp86_sto3g",
  paths_with_names=
      [
        (Path("C:/tmp/Molli/data/mol24_ex8_gfn2_crest_best_BP86_STO3G_opt_opt_steps.xyz"), "mol24_ex8_gfn2_crest_best_BP86_STO3G_opt_opt_steps"),
      ]
)


#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------
