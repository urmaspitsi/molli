# %%
trajectories_dir = Path("C:/tmp/gaussian/optimization_trajectories/mol24_ex16_crest5")
num_steps = 5

base_trajectory = MultiItemFileSource(
    file_path=trajectories_dir.joinpath("mol24_ex16_crest5_PBE1PBE_cc_pVTZ_TZVPFit_opt_steps.xyz"),
    name="pbe1pbe_cc_pvtz_tzvpfit",
    item_idxs=num_steps
  )

trajectories_to_compare = [
    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex16_gfn2_bp86_sto3g_crest_5_opt_steps.xyz"),
      name="bp86_sto3g",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex16_gfn2_bp86_def2svpp_svpfit_crest_5_opt_steps.xyz"),
      name="bp86_def2svpp_svpfit",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("crest_mol24_ex16_crest5_PBE1PBE_Def2SVPP_SVPFit_opt_steps.xyz"),
      name="pbe1pbe_def2svpp_svpfit",
      item_idxs=num_steps
      ),
  ]

metric_funcs = [
  (ms.mad_of_positions, "mean_abs_diff"),
  (ms.rmsd_of_positions, "root_mean_of_sq_diffs"),
]

mol24_ex16_crest5_traj_analyzer = TrajectoryAnalyzer(
    description="mol24_ex16_crest5 optimization trajectories",
    base_trajectory=base_trajectory,
    trajectories_to_compare=trajectories_to_compare,
    metric_functions=metric_funcs,
  )


# %%
trajectories_dir = Path("C:/tmp/gaussian/optimization_trajectories/mol24_ex16_crest5")
num_steps = 5

base_trajectory = MultiItemFileSource(
    file_path=trajectories_dir.joinpath("mol24_ex16_crest5_PBE1PBE_cc_pVTZ_TZVPFit_opt_steps.xyz"),
    name="pbe1pbe_cc_pvtz_tzvpfit",
    item_idxs=num_steps
  )
#%%
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import ase_utils as au
import gaussian_utils as GU
import metrics as ms
import utils as ut

from trajectory_analyzer import TrajectoryAnalyzer
from file_source import MultiItemFileSource
import metrics as ms


metric_funcs = [
  (ms.mad_of_positions, "mean_abs_diff"),
  (ms.rmsd_of_positions, "root_mean_of_sq_diffs"),
]

#%%
###############################################################################
# mol24_ex16_crest5
###############################################################################

trajectories_dir = Path("C:/tmp/gaussian/optimization_trajectories/mol24_ex16_crest5")

base_trajectory = MultiItemFileSource(
    file_path=trajectories_dir.joinpath("mol24_ex16_crest5_PBE1PBE_cc_pVTZ_TZVPFit_opt_steps.xyz"),
    name="pbe1pbe_cc_pvtz_tzvpfit",
    item_idxs=num_steps
  )

trajectories_to_compare = [
    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex16_gfn2_bp86_sto3g_crest_5_opt_steps.xyz"),
      name="bp86_sto3g",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex16_gfn2_bp86_def2svpp_svpfit_crest_5_opt_steps.xyz"),
      name="bp86_def2svpp_svpfit",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("crest_mol24_ex16_crest5_PBE1PBE_Def2SVPP_SVPFit_opt_steps.xyz"),
      name="pbe1pbe_def2svpp_svpfit",
      item_idxs=num_steps
      ),
  ]

mol24_ex16_crest5_traj_analyzer = TrajectoryAnalyzer(
    description="mol24_ex16_crest5 optimization trajectories",
    base_trajectory=base_trajectory,
    trajectories_to_compare=trajectories_to_compare,
    metric_functions=metric_funcs,
  )

# %%
mol24_ex16_crest5_traj_analyzer.metrics["mean_abs_diff"]["base_to_start"]


# %%
###############################################################################
# mol24_ex19_crest23
###############################################################################

trajectories_dir = Path("C:/tmp/gaussian/optimization_trajectories/mol24_ex19_crest23")

base_trajectory = MultiItemFileSource(
    file_path=trajectories_dir.joinpath("mol24_ex19_crest23_pbe1pbe_cc_pvtz_tzvpfit_opt_steps.xyz"),
    name="pbe1pbe_cc_pvtz_tzvpfit",
    item_idxs=num_steps
  )

trajectories_to_compare = [
    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex19_gfn2_bp86_sto3g_crest_23_opt_steps.xyz"),
      name="bp86_sto3g",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("crest_mol24_ex19_crest23_bp86_def2svpp_svpfit_opt_steps.xyz"),
      name="bp86_def2svpp_svpfit",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("crest_mol24_ex19_crest23_pbe1pbe_def2svpp_svpfit_opt_steps.xyz"),
      name="pbe1pbe_def2svpp_svpfit",
      item_idxs=num_steps
      ),
  ]


mol24_ex19_crest23_traj_analyzer = TrajectoryAnalyzer(
    description="mol24_ex19_crest23 optimization trajectories",
    base_trajectory=base_trajectory,
    trajectories_to_compare=trajectories_to_compare,
    metric_functions=metric_funcs,
  )

#%%
###############################################################################
# mol24_ex0a_c10
###############################################################################

trajectories_dir = Path("C:/tmp/gaussian/optimization_trajectories/mol24_ex0a_c10")

base_trajectory = MultiItemFileSource(
    file_path=trajectories_dir.joinpath("mol24_ex0a_c10_PBE1PBE_cc_pVTZ_TZVPFit_opt_steps.xyz"),
    name="pbe1pbe_cc_pvtz_tzvpfit",
    item_idxs=num_steps
  )

trajectories_to_compare = [
    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex0a_c10_bp86def2svpp_step10_PBE1PBE_cc_pVTZ_TZVPFit_opt_steps.xyz"),
      name="pbe1pbe_cc_pvtz_tzvpfit_preopt_bp86def2svpp_10steps",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex0a_gfn2_bp86_sto3g_crest_10_opt_steps.xyz"),
      name="bp86_sto3g",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex0a_gfn2_bp86_def2svpp_svpfit_crest_10_opt_steps.xyz"),
      name="bp86_def2svpp_svpfit",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex0a_gfn2_pbe1pbe_sto3g_crest_10_opt_steps.xyz"),
      name="pbe1pbe_sto3g",
      item_idxs=num_steps
      ),

    MultiItemFileSource(
      file_path=trajectories_dir.joinpath("mol24_ex0a_gfn2_pbe1pbe_def2svpp_svpfit_crest_10_opt_steps.xyz"),
      name="pbe1pbe_def2svpp_svpfit",
      item_idxs=num_steps
      ),
  ]

mol24_ex0a_crest10_traj_analyzer = TrajectoryAnalyzer(
    description="mol24_ex0a_crest10 optimization trajectories",
    base_trajectory=base_trajectory,
    trajectories_to_compare=trajectories_to_compare,
    metric_functions=metric_funcs,
  )

#%%
###############################################################################
#
###############################################################################

#%%

#%%

#%%


