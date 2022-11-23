import copy
from datetime import datetime
from itertools import chain
import json
import numpy as np
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

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
    new_lines = [xyz_parser.convert_xyz_coords_to_str(
        element=line[0],
        x=float(line[2].strip()),
        y=float(line[3].strip()),
        z=float(line[4].strip())
      ) for line in new_lines]

    num_atoms = len(new_lines)
    description = f"final xyz from {file_path.name}"

    res[results_key] = [f"{num_atoms}"] + [description] + new_lines

  except Exception as ex:
    res[error_key] = str(ex)

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
    description = f"opt step {i + 1} of {num_steps}"
    xyz = extract_opt_step_as_xyz_lines(
                                        lines=lines,
                                        between_lines=(start_line, end_line),
                                        description=description
                                        )
    res.extend(xyz) if collect_to_single_list else res.append(xyz)

  return res


def extract_scf_summary(
                        file_path: Union[str, Path],
                        collect_to_single_list: bool
                        ) -> List[str]:

  result_summary_dict = {}
  lines = ut.read_text_file_as_lines(file_path=file_path)
  block_lines = ut.get_block_start_line_nrs(
                                            lines=lines,
                                            search_text="SCF Done:"
                                            )

  num_steps = len(block_lines)
  text_lines = []
  separator_row = "-" * 90

  try:
    energy_start = float(lines[block_lines[0]].split("  ")[2])
    energy_end = float(lines[block_lines[-1]].split("  ")[2])
    energy_delta = str(round(energy_end - energy_start, 8)).rstrip()
  except Exception as ex:
    energy_delta = str(ex)

  try:
    minutes_per_step = 0
    total_minutes = 0
    elapsed_time_line_nrs = ut.get_block_start_line_nrs(lines, "Elapsed time:")
    if elapsed_time_line_nrs == None or len(elapsed_time_line_nrs) < 1:
      elapsed_time = "Elapsed time info not found"
      result_summary_dict["error"] = "Did not find the phrase: 'Elapsed time:' in the input file."
    else:
      time_line_idx = elapsed_time_line_nrs[0]
      elapsed_time = lines[time_line_idx].lstrip().rstrip().replace("time:      ", "time:").replace(" days", "d").replace(" hours", "h").replace(" minutes", "m").replace(" seconds.", "s")
      _,_, days, hours, minutes, seconds = elapsed_time.split()

      total_minutes = int(days[:-1].strip()) * 24 * 60 \
                    + int(hours[:-1].strip()) * 60 \
                    + int(minutes[:-1].strip()) \
                    + float(seconds[:-1].strip()) / 60
      minutes_per_step = round(total_minutes / num_steps, 1)

    result_summary_dict["elapsed_time_str"] = elapsed_time
    result_summary_dict["elapsed_time_minutes"] = round(total_minutes, 1)
    result_summary_dict["num_steps"] = num_steps
    result_summary_dict["minutes_per_step"] = minutes_per_step
    result_summary_dict["energy_start"] = energy_start
    result_summary_dict["energy_end"] = energy_end
    result_summary_dict["energy_delta"] = energy_delta

  except Exception as ex:
    elapsed_time = str(ex)
    result_summary_dict["error"] = str(ex)

  scf_data = [separator_row, f"SCF: change in energy = {energy_delta}   {elapsed_time}   {minutes_per_step} min/step", separator_row] \
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


def write_scf_summary_from_gaussian_logfile(
                                    log_file_path: Union[str, Path],
                                    output_path: Union[str, Path]
                                    ) -> Dict:

  try:
    summary = extract_scf_summary(
                                  file_path=log_file_path,
                                  collect_to_single_list=True
                                  )

    write_result = ut.write_text_file_from_lines(
                                  file_path=output_path,
                                  lines=summary["text_lines"]
                                )

    res = summary["summary"]
    res["write_result"] = write_result
    return res

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
                    ) -> Dict:

  out_dir = Path(input_path).parent \
              if output_dir == None \
              else Path(output_dir)

  input_file_name_stem = Path(input_path).stem

  output_file_name_stem = input_file_name_stem \
                            if len(input_file_name_stem) > 3 \
                            else Path(input_path).parent.parts[-1]

  task_results = {}

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

  output_path_scf = Path(out_dir).joinpath(f"{output_file_name_stem}_scf_summary.txt")
  task_results["scf_summary"] = write_scf_summary_from_gaussian_logfile(
                                    log_file_path=input_path,
                                    output_path=output_path_scf
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
                            aggregate_log_file_name: str="aggregate_log.txt"
                            ):

  res = [process_one_log_file(x, output_dir) for x in input_paths]

  # try sort ascending by final energy
  try:
    res.sort(key= lambda x: x["results"]["scf_summary"]["energy_end"])
  except:
    pass


  if output_dir != None:
    ut.write_text_file_json(
        file_name=Path(output_dir).joinpath(aggregate_log_file_name),
        data=res
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
