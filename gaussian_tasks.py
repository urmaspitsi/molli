
#%%

from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import gaussian_utils as GU
import utils as ut


#%%

############################################################################
# Create Gaussian job files: convert from .xyz trajectory (opt steps) file.
############################################################################

#xyz_file = Path("C:/tmp/gaussian/crest_mol24_ex15_BP86_STO3G/mol24_ex15_gfnff_crest_conformers.xyz")
#xyz_file = Path("C:/tmp/gaussian/crest_mol24_ex15_BP86_STO3G/mol24_ex16_gfn2_crest_conformers.xyz")

#xyz_file = Path("C:/tmp/gaussian/input/mol24_ex16_gfn2_crest_conformers.xyz")
#xyz_file = Path("C:/tmp/gaussian/input/mol24_ex15_gfnff_crest_conformers.xyz")
xyz_file = Path("C:/tmp/gaussian/input/mol24_ex19_gfn2_crestconfs_from_ex16_crest5.xyz")

gaussian_jobs_output_dir = Path("C:/tmp/gaussian/test_jobs") #

lines_to_prepend = [
    "#T BP86/Def2SVPP/SVPFit opt=(calcfc,maxcycles=5) formcheck",
    "",
    "mol24_ex19_gfn2_crestconfs_from_ex16_crest5",
    "",
    "0 1",
  ]

gaussian_job_files = GU.create_gaussian_job_files_from_xyz_steps(
                            input_path=xyz_file,
                            output_dir=gaussian_jobs_output_dir,
                            lines_before_xyz_coords=lines_to_prepend,
                            step_nrs_to_write=list(range(1, 62)),
                            job_file_name_prefix="mol24_ex19_gfn2_bp86_def2svpp_svpfit_crest_"
                          )


#%%
############################################################################
# Process Gaussian log files: log files as a list of full paths.
# Case 1: Many log files in one folder
############################################################################

log_files_dirs = [
#    Path("C:/tmp/gaussian/crest_reoptimize/crest_mol24_ex15_BP86_STO3G_results"),
#    Path("C:/tmp/gaussian/crest_reoptimize/crest_mol24_ex15_BP86_STO3G_results"),
    Path("C:/tmp/gaussian/crest_reoptimize/crest_mol24_ex15_BP86_Def2SVPP_SVPFit_results"),
  ]

gaussian_log_files = ut.get_file_paths_in_many_dirs(log_files_dirs, ".log")

#output_dir = "C:/tmp/gaussian/crest_reoptimize/logs"
output_dir = "C:/tmp/gaussian/crest_reoptimize/crest_mol24_ex15_BP86_Def2SVPP_SVPFit_logs"

#aggregate_log_file_name = "aggregate.log"
#aggregate_log_file_name = "ex15_bp86_sto3g_step10.log"
aggregate_log_file_name = "ex15_bp86_def2svpp_svpfit.log"

GU.process_many_log_files(
                          input_paths=gaussian_log_files,
                          output_dir=output_dir,
                          aggregate_log_file_name=aggregate_log_file_name,
#                          extract_summary_step_nr=10,
#                          do_only_summary=True
                          )


#%%
############################################################################
# Process Gaussian log files: multiple source folders, same log file name.
# Case 2: Many log files in different folders but log file name is same.
############################################################################

dir_names = [
    "mol5_B3P86_Def2SVPP_opt_24cpu",
    "mol5_PBE1PBE_SVP_SVPFit_opt_48cpu",
    "mol5_BP86_Def2SVPP_SVPFit_opt_24cpu",
    "mol24_BP86_STO3G_opt_48cpu",
    "mol24_BP86_STO3G_opt_fixOTi",
    "mol24_ex8_gfn2_crest_best_BP86_STO3G_opt",
    "mol24_ex8_gfn2_crest_best_PBE1PBE_STO3G_opt",
    "mol24_PBE1PBE_STO3G_opt",
  ]


gaussian_log_files = [Path(f"C:/tmp/gaussian/jobs_results/{x}/job.log") for x in dir_names]

output_dir = "C:/tmp/gaussian/jobs_logs"

aggregate_log_file_name = "aggregate.log"

GU.process_many_log_files(
                          input_paths=gaussian_log_files,
                          output_dir=output_dir,
                          aggregate_log_file_name=aggregate_log_file_name
                          )


#%%
############################################################################
# Process specific Gaussian log files
############################################################################

gaussian_log_files = [
  Path("C:/tmp/gaussian/jobs_results/mol24_PBE1PBE_cc_pVTZTZVPFit_energy/job.log"),
  Path("C:/tmp/gaussian/jobs_results/mol24_PBE1PBE_cc_pVTZTZVPFit_energy_g02/job.log"),
  Path("C:/tmp/gaussian/jobs_juhan/mol5_mol24_PBE1PBE_cc_pVTZ_TZVPFit/mol5_PBE1PBE_cc_pVTZ_TZVPFit_job.log"),
  Path("C:/tmp/gaussian/jobs_juhan/mol5_mol24_PBE1PBE_cc_pVTZ_TZVPFit/mol24_PBE1PBE_cc_pVTZ_TZVPFit_job.log"),
  Path("C:/tmp/gaussian/jobs_juhan/mol5_mol24_PBE1PBE_cc_pVTZ_TZVPFit/mol5_PBE1PBE_cc_pVTZ_TZVPFit_opt_job.log"),
  Path("C:/tmp/gaussian/jobs_juhan/mol5_mol24_PBE1PBE_cc_pVTZ_TZVPFit/mol24_PBE1PBE_cc_pVTZ_TZVPFit_opt_job.log"),

  ]

output_dir = "C:/tmp/gaussian/jobs_logs2"

aggregate_log_file_name = "aggregate.log"

GU.process_many_log_files(
                          input_paths=gaussian_log_files,
                          output_dir=output_dir,
                          aggregate_log_file_name=aggregate_log_file_name
                          )

#%%
# Test extracting final xyz

gaussian_log_file = Path("C:/tmp/gaussian/jobs_results/mol24_PBE1PBE_cc_pVTZTZVPFit_energy/job.log")

output_path = "C:/tmp/gaussian/jobs_logs2/aaa_xyz_test.txt"

xyz_lines = GU.extract_final_xyz(file_path=gaussian_log_file)["results"]

ut.write_text_file_from_lines(
                              file_path=output_path,
                              lines=xyz_lines
                            )


# %%
len(xyz_lines)
# %%
xyz_lines
# %%
