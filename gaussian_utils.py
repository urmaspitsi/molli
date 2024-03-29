import copy
from datetime import datetime
from itertools import chain
import json
import numpy as np
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import ase_utils as au
import constants as C
import utils as ut
import xyz_parser

################################################################
# Gaussian logfile parsing
################################################################


def write_xyz_from_gaussian_logfile(
                                    log_file_path: Union[str, Path],
                                    output_path: Union[str, Path]
                                    ) -> int:

  xyz_from_gaussian_log = extract_final_xyz(log_file_path)

  results_key = "results"
  error_key = "error"

  is_error = True if error_key in xyz_from_gaussian_log \
                    and xyz_from_gaussian_log[error_key] != None \
                    and len(xyz_from_gaussian_log[error_key]) > 0 \
                  else False

  has_results = True if results_key in xyz_from_gaussian_log \
                        and xyz_from_gaussian_log[results_key] != None \
                        and len(xyz_from_gaussian_log[results_key]) > 1 \
                     else False

  if has_results:
    return ut.write_text_file_from_lines(
      file_path=output_path,
      lines=xyz_from_gaussian_log[results_key]
    )

  elif is_error:
    return xyz_from_gaussian_log[error_key]
  else:
    return xyz_from_gaussian_log


def get_start_end_line_nr_for_xyz_block(
                                        lines: List[str],
                                        between_lines: Tuple[int,int]
                                        ) -> Tuple[int,int]:

  ''' xyz block
      return (start_line_nr, end_line_nr)
  '''
  start_nr, end_nr = between_lines
  start_nr +=5 # this should be the first line of the xyz block
  search_txt = "--------------------------------------------"

  for i in range(start_nr, end_nr):
    if search_txt in lines[i]:
      break

  return (start_nr, i)


def prepend_xyz_info(block: List[str], description: str) -> List[str]:
  return [str(len(block))] + [description] + block


def transform_xyz_block_row(row: str) -> str:
  _, atomic_nr, _, x, y, z = row.lstrip().rstrip().split()[:6]
  element = C.atomic_numbers_to_elements[int(atomic_nr)]

  len_el = len(element)
  len_x = len(x)
  len_y = len(y)
  len_z = len(z)

  space_1 = " " * (20 - len_el - len_x)
  space_2 = " " * (38 - len_el - len(space_1) - len_x - len_y)
  space_3 = " " * (56 - len_el - len(space_1) - len_x - len(space_2) - len_y - len_z)

  return f"{element}{space_1}{x}{space_2}{y}{space_3}{z}"


def extract_opt_step_as_xyz_lines(
                                    lines: List[str],
                                    between_lines: Tuple[int,int],
                                    description: str
                                  ) -> List[str]:

  block_bounds = get_start_end_line_nr_for_xyz_block(
                                    lines=lines,
                                    between_lines=between_lines
                                  )

  block = [lines[i] for i in range(block_bounds[0], block_bounds[1])]

  transformed_block = [transform_xyz_block_row(i) for i in block]

  final_xyz_block = prepend_xyz_info(
                                      transformed_block,
                                      description=description
                                    )

  return final_xyz_block


def extract_elapsed_time(lines: List[str]) -> str:
  res = ""

  try:
    line_nr = ut.get_block_start_line_nrs(
                                          lines=lines,
                                          search_text="Elapsed time:"
                                          )[0]

    res = lines[line_nr].strip().replace(" days", "d").replace(" hours", "h").replace(" minutes", "m").replace(" seconds.", "s")
    res = res.split(":")[1].strip()
  except:
    res = ""

  return res


def extract_job_completion_datetime(lines: List[str]) -> str:
  res = ""

  try:
    # completion time usually last line
    res = lines[-1].strip()
    if not "Normal termination of Gaussian" in res:
      # try error termination
      line_nr = ut.get_block_start_line_nrs(
                                            lines=lines,
                                            search_text="Error termination via"
                                            )[0]
      res = lines[line_nr]

    res = res.split(" at ")[1]
    res = res.replace(".", "").strip()
    weekday, month_name, day, time, year = res.split()
    res = f"{day}-{month_name}-{year} {time}"

  except Exception as ex:
    res = str(ex)

  return res


def extract_job_cpu_time(lines: List[str]) -> str:
  res = ""

  try:
    line_nr = ut.get_block_start_line_nrs(
                                          lines=lines,
                                          search_text="Job cpu time:"
                                          )[0]

    res = lines[line_nr].strip().replace(" days", "d").replace(" hours", "h").replace(" minutes", "m").replace(" seconds.", "s")
    res = res.split(":")[1].strip()

  except:
    res = ""

  return res


def get_total_minutes_from_elapsed_time(elapsed_time: str) -> int:
  if elapsed_time == None or len(elapsed_time) < 1:
    return 0
  else:
    total_minutes = 0

    try:
      days, hours, minutes, seconds = elapsed_time.split()

      total_minutes = int(days[:-1].strip()) * 24 * 60 \
                    + int(hours[:-1].strip()) * 60 \
                    + int(minutes[:-1].strip()) \
                    + float(seconds[:-1].strip()) / 60

      total_minutes = round(total_minutes, 1)
    except:
      total_minutes = 0

    return total_minutes


def extract_final_xyz(
                      file_path: Union[str, Path],
                      ) -> List[str]:

  results_key = "results"
  error_key = "error"
  res = {
      results_key: "",
      error_key: "",
      "results_key": results_key,
      "error_key": error_key
    }

  try:

    lines = ut.read_text_file_as_lines(file_path=file_path)

    start_line = ut.get_block_start_line_nrs(
                  lines=lines,
                  search_text="Unable to Open any file for archive entry."
                  )[0]

    end_line = ut.get_block_start_line_nrs(
                  lines=lines,
                  search_text="The archive entry for this job was punched."
                  )[0]

    new_str = "".join([lines[i].strip() for i in range(start_line, end_line)])
    new_lines = new_str.split("\\")
    new_lines = [x.strip() for x in new_lines]
    
    xyz_start_line = ut.get_line_nrs_starts_with_text(
                                              lines=new_lines,
                                              search_text="0,1"
                                              )[0]

    xyz_end_line = ut.get_line_nrs_starts_with_text(
                                            lines=new_lines,
                                            search_text="Version="
                                            )[0]

    new_lines = new_lines[(xyz_start_line + 1):xyz_end_line]
    new_lines = [line.split(",") for line in new_lines if len(line) > 8]
    skip_2nd_col = True if len(new_lines[0]) > 4 else False # energy logs have atomic nr in second col.
    col_x = 2 if skip_2nd_col else 1
    new_lines = [xyz_parser.convert_xyz_coords_to_str(
        element=line[0],
        x=float(line[col_x].strip()),
        y=float(line[col_x + 1].strip()),
        z=float(line[col_x + 2].strip())
      ) for line in new_lines]

    num_atoms = len(new_lines)
    description = f"final xyz from {file_path.name}"

    res[results_key] = [f"{num_atoms}"] + [description] + new_lines

  except Exception as ex:
    res[error_key] = str(ex)

  return res


def extract_energy(lines: List[str]) -> str:
  res = ""

  try:
    res = lines[ut.get_block_start_line_nrs(
                                            lines=lines,
                                            search_text="SCF Done:"
                                            )[0]]

    res = res.split("=")[1].strip()
    res = res.split()[0].strip()

  except:
    res = ""

  return res


def extract_optimization_steps_as_xyz(
                                      file_path: Union[str, Path],
                                      collect_to_single_list: bool
                                      ) -> List[str]:

  lines = ut.read_text_file_as_lines(file_path=file_path)
  xyz_block_lines = ut.get_block_start_line_nrs(
                                                lines=lines,
                                                search_text="Standard orientation:"
                                                )

  num_steps = len(xyz_block_lines)
  res = []
  for i in range(num_steps):
    start_line = xyz_block_lines[i]
    end_line = xyz_block_lines[i + 1] if i < num_steps - 1 else len(lines)
    energy = extract_energy(lines=lines[start_line:end_line])
    if len(energy) < 1 and i > 0:
      # try previous: if converged then last two xyz blocks are identical and SCF Done is not shown for the last.
      energy = extract_energy(lines=lines[xyz_block_lines[i - 1]:start_line])

    description = f"{energy}, opt step {i + 1} of {num_steps}, source: {file_path.name}"
    xyz = extract_opt_step_as_xyz_lines(
                                        lines=lines,
                                        between_lines=(start_line, end_line),
                                        description=description
                                        )
    res.extend(xyz) if collect_to_single_list else res.append(xyz)

  return res


def extract_chemical_formula(lines: List[str]) -> str:
  res = ""

  try:
    res = lines[ut.get_block_start_line_nrs(
                                            lines=lines,
                                            search_text="Framework group"
                                            )[0]].strip()

    res = res.split("(")[1]
    res = res.split(")")[0].strip()

  except:
    res = ""

  return res


def extract_gaussian_version(lines: List[str]) -> str:
  res = ""

  try:
    res = lines[ut.get_block_start_line_nrs(
                                            lines=lines,
                                            search_text="Cite this work as:"
                                            )[0] + 1].strip()

  except:
    res = ""

  return res


def extract_dft_info(lines: List[str]) -> str:
  res = ""

  try:
    dft_info_line_nr = ut.get_block_start_line_nrs(
                                              lines=lines,
                                              search_text="Rotational constants"
                                              )[0] + 1
    basis_functions_info = lines[dft_info_line_nr].strip()
    electrons_info = lines[dft_info_line_nr + 1].strip()
    res = f"{basis_functions_info}, {electrons_info}"

  except:
    res = ""

  return res


def extract_gaussian_command_and_dft_functional(lines: List[str]) -> Tuple[str, str]:
  res = ("", "")

  try:
    gaussian_command = lines[ut.get_block_start_line_nrs(
                                              lines=lines,
                                              search_text=" #"
                                              )[0]].strip()
    dft_functional = gaussian_command.split()[1]
    res = (gaussian_command, dft_functional)

  except:
    res = ("", "")

  return res


def extract_scf_summary(
                        file_path: Union[str, Path],
                        collect_to_single_list: bool,
                        max_step_nr: int=0,
                        ignore_shorter_runs: bool=False,
                        return_converged_only: bool=False
                        ) -> List[str]:

  skip = False # if error is found then skip subsequent code block

  result_summary_dict = {}
  lines = ut.read_text_file_as_lines(file_path=file_path)

  result_summary_dict["gaussian_version"] = extract_gaussian_version(lines=lines)

  gaussian_command, dft_functional = extract_gaussian_command_and_dft_functional(lines=lines)
  result_summary_dict["gaussian_command"] = gaussian_command
  result_summary_dict["dft_functional"] = dft_functional

  result_summary_dict["dft_info"] = extract_dft_info(lines=lines)
  result_summary_dict["chemical_formula"] = extract_chemical_formula(lines=lines)

  block_lines = ut.get_block_start_line_nrs(
                                            lines=lines,
                                            search_text="SCF Done:"
                                            )

  converged_line_nrs = ut.get_block_start_line_nrs(
                                            lines=lines,
                                            search_text="Stationary point found"
                                            )

  is_converged = True \
    if converged_line_nrs and len(converged_line_nrs) > 0 \
    else False

  result_summary_dict["optimization_converged"] = is_converged

  if return_converged_only and not is_converged:
      result_summary_dict["error"] = "Optimization not converged"
      skip = True

  # Collect summary at the max_step_nr
  num_steps_original = len(block_lines)
  if len(block_lines) > 1 and max_step_nr > 0:
    block_lines = block_lines[:max_step_nr]

    if ignore_shorter_runs and len(block_lines) < max_step_nr:
      result_summary_dict["error"] = f"Too few optimization steps, num_steps_required={max_step_nr}, num_steps_available={len(block_lines)}"
      skip = True

  num_steps = len(block_lines)
  steps_pct = 1 if num_steps_original == num_steps else num_steps / num_steps_original

  text_lines = []
  separator_row = "-" * 90
  energy_delta = 9999 # some arbitrary large number, if error

  if not skip:

    try:
      energy_start = float(lines[block_lines[0]].split("  ")[2])
      energy_end = float(lines[block_lines[-1]].split("  ")[2])
      energy_delta = energy_end - energy_start
      energy_delta_str = f"{round(energy_delta, 8)} a.u., {round(energy_delta * C.hartree_in_kcal_per_mol, 2)} kcal/mol, {round(energy_delta * C.hartree_in_kJ_per_mol, 2)} kJ/mol"

    except Exception as ex:
      energy_delta_str = str(ex)

    try:
      minutes_per_step = 0
      total_minutes = 0
      elapsed_time = extract_elapsed_time(lines=lines)
      if elapsed_time == None or len(elapsed_time) < 1:
        elapsed_time = "Elapsed time info not found"
        result_summary_dict["error"] = "Did not find the phrase: 'Elapsed time:' in the input file."
      else:
        total_minutes = get_total_minutes_from_elapsed_time(elapsed_time=elapsed_time)
        total_minutes *= steps_pct # scale total minutes, if partial summary
        minutes_per_step = round(total_minutes / num_steps, 1)

      elapsed_time_str = elapsed_time \
                          if steps_pct == 1 \
                          else f"{elapsed_time} (time for {num_steps_original} steps)"

      result_summary_dict["elapsed_time_str"] = elapsed_time_str
      result_summary_dict["elapsed_time_minutes"] = round(total_minutes, 1)
      result_summary_dict["num_steps"] = num_steps
      result_summary_dict["minutes_per_step"] = minutes_per_step
      result_summary_dict["energy_start"] = energy_start
      result_summary_dict["energy_end"] = energy_end
      result_summary_dict["energy_delta"] = energy_delta
      result_summary_dict["energy_delta_text"] = energy_delta_str

      job_cpu_time = extract_job_cpu_time(lines=lines)
      job_cpu_minutes = get_total_minutes_from_elapsed_time(elapsed_time=job_cpu_time)

      result_summary_dict["job_cpu_time"] = job_cpu_time
      result_summary_dict["job_cpu_hours"] = round(job_cpu_minutes / 60, 1)
      result_summary_dict["job_cpu_hours_per_step"] = round(job_cpu_minutes / num_steps / 60, 2)
      result_summary_dict["job_completion_datetime"] = extract_job_completion_datetime(lines=lines)

    except Exception as ex:
      elapsed_time = str(ex)
      result_summary_dict["error"] = str(ex)

    scf_data = [separator_row, f"SCF: change in energy = {energy_delta_str}  {elapsed_time_str}  {minutes_per_step} min/step", separator_row] \
      + [f"{i + 1}:" + " " * (4 - len(str(i + 1))) + lines[x].replace("\n", "") for i,x in enumerate(block_lines)] \
      + ["\n"]

    text_lines.extend(scf_data) if collect_to_single_list else text_lines.append(scf_data)

    for i in range(num_steps):
      start_line = block_lines[i]
      for j in range(100):
        start_line -= 1
        if "------------------------" in lines[start_line] or start_line < 0:
          start_line += 1
          break

      end_line = block_lines[i] + 3
      description = f"opt step {i + 1} of {num_steps}"
      summary = [separator_row, description, separator_row] + lines[start_line:end_line] + ["\n"]
      summary = [x.replace("\n", "") for x in summary]
      text_lines.extend(summary) if collect_to_single_list else text_lines.append(summary)

    # final message: wall time etc
    final_block = [separator_row, "END", separator_row] + lines[-13:]
    final_block = [x.replace("\n", "") for x in final_block]
    text_lines.extend(final_block) if collect_to_single_list else text_lines.append(final_block)

  res = {}
  res["summary"] = result_summary_dict
  res["text_lines"] = text_lines
  return res


def extract_and_write_scf_summary_from_gaussian_logfile(
                                    log_file_path: Union[str, Path],
                                    output_path: Union[str, Path, None],
                                    max_step_nr: int=0,
                                    ignore_shorter_runs: bool=False,
                                    return_converged_only: bool=False
                                    ) -> Dict:

  try:
    summary = extract_scf_summary(
                                  file_path=log_file_path,
                                  collect_to_single_list=True,
                                  max_step_nr=max_step_nr,
                                  ignore_shorter_runs=ignore_shorter_runs,
                                  return_converged_only=return_converged_only
                                  )

    first_part = {
      "scf_summary_file": "",
      }

    if output_path:
      write_result = ut.write_text_file_from_lines(
                                    file_path=output_path,
                                    lines=summary["text_lines"]
                                  )

      first_part = {
        "scf_summary_file": Path(output_path).name,
        }

    last_part = summary["summary"]
    res = {**first_part, **last_part}
    return res

  except Exception as ex:
    return str(ex)


def write_optimization_steps_from_gaussian_logfile(
                                    log_file_path: Union[str, Path],
                                    output_path: Union[str, Path]
                                    ) -> int:
  try:
    xyz_blocks = extract_optimization_steps_as_xyz(
                                                  file_path=log_file_path,
                                                  collect_to_single_list=True
                                                  )

    return ut.write_text_file_from_lines(file_path=output_path, lines=xyz_blocks)

  except Exception as ex:
    return str(ex)


def write_last_optimization_step_from_gaussian_logfile(
                                    log_file_path: Union[str, Path],
                                    output_path: Union[str, Path]
                                    ) -> int:
  try:
    xyz_blocks = extract_optimization_steps_as_xyz(
                                                  file_path=log_file_path,
                                                  collect_to_single_list=False
                                                  )

    last_block = xyz_blocks[-1]
    return ut.write_text_file_from_lines(file_path=output_path, lines=last_block)

  except Exception as ex:
    return str(ex)


def get_gaussian_log_file_type(file_path: Union[str, Path]) -> str:

  lines = ut.read_text_file_as_lines(file_path=file_path)

  optimization_line_nrs = ut.get_block_start_line_nrs(
        lines=lines,
        search_text="GradGradGradGradGradGradGradGradGradGradGrad"
      )

  res = "energy"
  if optimization_line_nrs and len(optimization_line_nrs) > 0:
    res = "optimization"

  return res


def process_one_log_file(
                    input_path: Union[str, Path],
                    output_dir: Union[str, Path]=None,
                    extract_summary_step_nr: int=0,
                    ignore_shorter_runs: bool=False,
                    do_only_summary: bool=False,
                    return_converged_only: bool=False
                    ) -> Dict:

  out_dir = Path(input_path).parent \
              if output_dir == None \
              else Path(output_dir)

  input_file_name_stem = Path(input_path).stem

  output_file_name_stem = input_file_name_stem \
                            if len(input_file_name_stem) > 3 \
                            else Path(input_path).parent.parts[-1]

  task_results = {}

  if not do_only_summary:
    gaussian_log_file_type = get_gaussian_log_file_type(file_path=input_path)

    output_path_xyz = Path(out_dir).joinpath(f"{output_file_name_stem}.xyz")
    task_results["final_xyz"] = write_xyz_from_gaussian_logfile(
                                      log_file_path=input_path,
                                      output_path=output_path_xyz
                                      )

    if gaussian_log_file_type.startswith("optim"):
      output_path_opt_steps = Path(out_dir).joinpath(f"{output_file_name_stem}_opt_steps.xyz")
      task_results["opt_steps"] = write_optimization_steps_from_gaussian_logfile(
                                        log_file_path=input_path,
                                        output_path=output_path_opt_steps
                                        )

      output_path_last_opt_step = Path(out_dir).joinpath(f"{output_file_name_stem}_last_step.xyz")
      task_results["last_opt_step"] = write_last_optimization_step_from_gaussian_logfile(
                                        log_file_path=input_path,
                                        output_path=output_path_last_opt_step
                                        )
    else:
      output_path_standard_orientation = Path(out_dir).joinpath(f"{output_file_name_stem}_standard_orientation.xyz")
      task_results["standard_orientation"] = write_last_optimization_step_from_gaussian_logfile(
                                        log_file_path=input_path,
                                        output_path=output_path_standard_orientation
                                        )

  output_path_scf = None \
                    if do_only_summary \
                    else Path(out_dir).joinpath(f"{output_file_name_stem}_scf_summary.txt")

  task_results["scf_summary"] = extract_and_write_scf_summary_from_gaussian_logfile(
                                    log_file_path=input_path,
                                    output_path=output_path_scf,
                                    max_step_nr=extract_summary_step_nr,
                                    ignore_shorter_runs=ignore_shorter_runs,
                                    return_converged_only=return_converged_only
                                    )

  res = {
      "input_path": str(input_path),
      "output_dir": str(out_dir),
      "results": task_results,
    }

  return res


def process_many_log_files(
                            input_paths: List[Union[str, Path]],
                            output_dir: Union[str, Path]=None,
                            aggregate_log_file_name: str="aggregate.log",
                            extract_summary_step_nr: int=0,
                            ignore_shorter_runs: bool=False,
                            do_only_summary: bool=False,
                            return_converged_only: bool=False,
                            write_last_opt_steps_file_path: Path=None,
                            show_errors_in_output: bool=True
                            ) -> Dict:

  '''
    Processes Gaussian log files.
    Returns Dictionary of with results.

    extract_summary_step_nr: step nr (or generally [item index + 1] xyz-item)
      from xyz file containing multiple xyz-data-sections, i.e. more than one geometries.
      default = 0, i.e. returns last xyz-item from the file.
      step_nr=1 corresponds to 0 in 0-based idx-s.

    ignore_shorter_runs: if Ture and step_nr > 1 then returns data that corresponds to step_nr.
     It might be useful if you want to compare trajectories that contain
     at least step_nr amount of steps.

    do_only_summary: if True then writes aggregate summary file.
                     if False then writes following information:
                        - aggregate summary file
                        - for each item in input_paths: 
                            - optimization steps as xyz file.
                            - last optimization step as xyz file.
                            - last xyz section from Gaussian log. It should be the same geometry as
                                  the last optimization step xyz file, just in different alignment.
                            - scf summary

      default = False (i.e. writes a lot of files)

    write_last_opt_steps_file_path: if specified then collects all xyz-items from the last-optimization step
      into one file, and writes it to the path specified (write_last_opt_steps_file_path).
      Might be useful for aggregating conformers from various optimization runs into one xyz file.
      default = None

    show_errors_in_output: whether occurred errors are returned in the final output.
      default = True

  '''

  res = [process_one_log_file(
            x,
            output_dir,
            extract_summary_step_nr=extract_summary_step_nr,
            ignore_shorter_runs=ignore_shorter_runs,
            do_only_summary=do_only_summary,
            return_converged_only=return_converged_only
          ) for x in input_paths]

  summary = {}
  errors = []

  # filter out experiments with missing data
  items_with_errors = []
  valid_items = []
  valid_paths = []
  for idx, res_item in enumerate(res):
    is_error = "error" in res_item["results"]["scf_summary"] \
        or "energy_end" not in res_item["results"]["scf_summary"]
    if is_error:
      res_item["aggregate_error"] = {}
      try:
        res_item["aggregate_error"]["input_path"] = res_item["input_path"]
        res_item["aggregate_error"]["error"] = res_item["results"]["scf_summary"]["error"]
        res_item["aggregate_error"]["scf_summary_file"] = res_item["results"]["scf_summary"]["scf_summary_file"]
        
      except:
        pass

      items_with_errors.append(res_item)
    else:
      valid_items.append(res_item)
      valid_paths.append(input_paths[idx])

  # items_with_errors = [x for x in res if "energy_end" not in x["results"]["scf_summary"]]
  if len(items_with_errors) > 0:
    errors.append({
      "missing_data_error": "Necessary data is missing from log file. Further details in the specific item's summary file.",
      "items_with_errors": [x["aggregate_error"] for x in items_with_errors],
      })
    res = valid_items

  # try sort ascending by final energy
  diff_best_worst_str = ""
  rank_list = []
  try:
    res.sort(key= lambda x: x["results"]["scf_summary"]["energy_end"])
    best_energy = res[0]["results"]["scf_summary"]["energy_end"]
    worst_energy = res[-1]["results"]["scf_summary"]["energy_end"]
    diff_best_worst = best_energy - worst_energy
    diff_best_worst_str = f"{round(diff_best_worst,6)} a.u., {round(diff_best_worst * C.hartree_in_kcal_per_mol, 2)} kcal/mol, {round(diff_best_worst * C.hartree_in_kJ_per_mol, 2)} kJ/mol"
    for i, dct in enumerate(res):
      energy_diff_to_best = dct["results"]["scf_summary"]["energy_end"] - best_energy
      inp_path = Path(dct["input_path"])
      file_name = inp_path.name
      name = file_name if len(file_name) > 10 else f"{inp_path.parent.name}_{file_name}"
      energy_diff_str = f"{round(energy_diff_to_best,6)} a.u., {round(energy_diff_to_best * C.hartree_in_kcal_per_mol, 2)} kcal/mol, {round(energy_diff_to_best * C.hartree_in_kJ_per_mol, 2)} kJ/mol"
      rank_list.append(f"{i + 1}: diff best: {energy_diff_str}, source: {name}")
      dct["results"]["scf_summary"]["energy_diff_to_best"] = energy_diff_str
  except Exception as ex:
    errors.append(str(ex))

  summary["num_experiments_total"] = len(res) + len(items_with_errors)
  summary["num_experiments_successful"] = len(res)
  summary["num_experiments_failed"] = len(items_with_errors)
  summary["energy_diff_best_worst"] = diff_best_worst_str
  summary["ranking"] = rank_list

  if write_last_opt_steps_file_path:
    valid_paths_sorted = {}
    try:
      valid_paths_sorted = {
        dct["input_path"]: dct["results"]["scf_summary"]["energy_end"] for dct in res
      }

    except:
      valid_paths_sorted = {}

    if len(valid_paths_sorted) > 0:
      valid_paths = [Path(k) for k, _ in valid_paths_sorted.items()]

    try:
      if extract_summary_step_nr > 0:
        last_opt_steps_all_runs = [extract_optimization_steps_as_xyz(
                                                      file_path=x,
                                                      collect_to_single_list=False
                                                      ) for x in valid_paths]

        # shorter runs should be already filtered out above and reflected in valid_paths
        # take last run or specified run if available
        last_opt_steps_xyz = [
                              x[-1] \
                                if len(x) < extract_summary_step_nr \
                                else x[extract_summary_step_nr - 1]
                              for x in last_opt_steps_all_runs
                            ]

      else:
        last_opt_steps_xyz = [extract_optimization_steps_as_xyz(
                                                      file_path=x,
                                                      collect_to_single_list=False
                                                      )[-1] for x in valid_paths]

      dicts = [xyz_parser.convert_xyz_lines_to_dict(x, convert_coords_to_float=True) for x in last_opt_steps_xyz]

      mols = [au.create_ase_atoms_from_xyz_data(xyz_data=xyz_data) for xyz_data in dicts]

      aggregate_xyz_res = au.write_ase_atoms_to_xyz_file(
                atoms_list=mols,
                output_path=write_last_opt_steps_file_path
              )

      aligned_mols = [
          au.align_2_molecules_min_rmsd(
                                      target=mols[0],
                                      atoms_to_align=x
                                    ) for x in mols
        ]

      out_file_stem = write_last_opt_steps_file_path.stem
      out_file_suffix = write_last_opt_steps_file_path.suffix
      aligned_file_name = f"{out_file_stem}_aligned{out_file_suffix}"
      aligned_mols_path = write_last_opt_steps_file_path.parent.joinpath(aligned_file_name)

      aggregate_xyz_res = au.write_ase_atoms_to_xyz_file(
                atoms_list=aligned_mols,
                output_path=aligned_mols_path
              )

      summary["last_opt_steps_file"] = aggregate_xyz_res

    except Exception as ex:
      errors.append(str(ex))

  if len(errors) > 0 and show_errors_in_output:
    summary["error"] = errors

  result_dict = {
                  "summary": summary,
                  "experiments": res,
                }

  if output_dir != None:
    ut.write_text_file_json(
        file_name=Path(output_dir).joinpath(aggregate_log_file_name),
        data=result_dict
      )


def create_gaussian_job_files_from_xyz_steps(
                              input_path: Union[str, Path],
                              output_dir: Union[str, Path],
                              lines_before_xyz_coords: List[str],
                              step_nrs_to_write: Union[None, List[int]],
                              job_file_name_prefix: Union[None, str],
                            ) -> List[Path]:

  xyz_steps = xyz_parser.read_xyz_file(
                                        input_path=input_path,
                                        convert_coords_to_float=False
                                      )

  res = []

  file_name_prefix = "job_" \
                        if job_file_name_prefix == None \
                          or len(job_file_name_prefix) < 1 \
                        else job_file_name_prefix

  for step_nr in step_nrs_to_write:
    if step_nr >= 0 and step_nr <= len(xyz_steps):
      job_file_name = Path(output_dir).joinpath(f"{file_name_prefix}{step_nr}.gjf")
      xyz = xyz_steps[step_nr - 1]
      name_idx = len(lines_before_xyz_coords) -  3

      lines_prepend = [f"{x}_{step_nr}" \
                         if i == name_idx \
                         else x for i,x in enumerate(lines_before_xyz_coords)]

      job_file_lines = lines_prepend + xyz["xyz_lines"] + [" "] # add space in the end just in case. Gaussian may crash if not sapce in the end ?!
      ut.write_text_file_from_lines(
                                  file_path=job_file_name,
                                  lines=job_file_lines
                                )
      res.append(job_file_name)

  return res











