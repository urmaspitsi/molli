
#%%
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import gaussian_utils as GU
import utils as ut


#%%
############################################################################
# Create Gaussian job files: convert from .xyz trajectory (opt steps) file.
############################################################################

xyz_file = Path("C:/tmp/gaussian/input/ex0ff_crestconfs.xyz")

gaussian_jobs_output_dir = Path("C:/tmp/gaussian/test_jobs") #

#    "#T PBE1PBE/Def2SVPP/SVPFit opt=(calcfc,maxcycles=10) formcheck",
#    "#T BP86/Def2SVPP/SVPFit opt=(calcfc,maxcycles=10) formcheck",
#    "#T BP86/Def2SVPP/SVPFit opt=(calcfc,maxcycles=100) formcheck",
#    "#T BP86/STO-3G opt=(calcfc,maxcycles=25) formcheck",
#    "#T PBE1PBE/STO-3G opt=(calcfc,maxcycles=25) formcheck",

lines_to_prepend = [
    "#T BP86/Def2SVPP/SVPFit opt=(calcfc,maxcycles=10) formcheck",
    "",
    "ex0ff_gfnff_crestconfs",
    "",
    "0 1",
  ]

gaussian_job_files = GU.create_gaussian_job_files_from_xyz_steps(
                            input_path=xyz_file,
                            output_dir=gaussian_jobs_output_dir,
                            lines_before_xyz_coords=lines_to_prepend,
                            step_nrs_to_write=list(range(1, 100)),
                            job_file_name_prefix="ex0ff_gfnff_bp86_def2svpp_svpfit_crest_"
                          )

# top5 from ex0a_BP86_STO3G: 13,29,23,6,26
# top5 from ex0a_BP86_Def2SVPP_SVPFit: 10,1,13,9,20
# combined top5: [1,6,9,10,13,20,23,26,29]

#%%
############################################################################
# 
############################################################################

# bp86_sto3g
# bp86_def2svpp_svpfit
# pbe1pbe_def2svpp_svpfit
# pbe1pbe_sto3g

#%%
############################################################################
# Process Gaussian log files: log files as a list of full paths.
############################################################################

# bp86_sto3g
# bp86_def2svpp_svpfit
# pbe1pbe_def2svpp_svpfit
# pbe1pbe_sto3g

log_files_dirs = [
    Path("C:/tmp/gaussian/optimize/crest_s4_bp86_sto3g/gaussian_results"),
  ]

gaussian_log_files = ut.get_file_paths_in_many_dirs(log_files_dirs, ".log")

output_dir = Path("C:/tmp/gaussian/optimize/crest_s4_bp86_sto3g/molli_results")

aggregate_log_file_name = "crest_s4_bp86_sto3g_log.txt"

aggregate_xyz_file = output_dir.joinpath("crest_s4_bp86_sto3g_confs.xyz")

GU.process_many_log_files(
                          input_paths=gaussian_log_files,
                          output_dir=output_dir,
                          aggregate_log_file_name=aggregate_log_file_name,
#                          extract_summary_step_nr=10,
#                          ignore_shorter_runs=False,
                          do_only_summary=True,
                          write_last_opt_steps_file_path=aggregate_xyz_file
                          )

#%%

