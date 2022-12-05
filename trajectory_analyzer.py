from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, NamedTuple, Set, Tuple, Union
from ase import Atoms

import ase_utils as au
from file_source import MultiItemFileSource

@dataclass
class TrajectoryAnalyzer():
  description: str
  base_trajectory: MultiItemFileSource
  trajectories_to_compare: List[MultiItemFileSource] = field(default_factory=list)
  metric_functions: List[Tuple[Callable, str]] = field(default_factory=list)

  base_at_steps: List[Atoms] = field(init=False)
  mols_to_compare_at_steps: List[List[Atoms]] = field(init=False)
  metrics: Dict = field(init=False)

  def __post_init__(self):
    self.base_at_steps = au.extract_mols_from_xyz_files(sources=[self.base_trajectory])[0]
    self.mols_to_compare_at_steps = au.extract_mols_from_xyz_files(sources=self.trajectories_to_compare)
    if self.metric_functions and len(self.metric_functions) > 0:
      self.metrics = self.calulate_all_metrics(metric_funcs=self.metric_functions)

  def __getitem__(self, i: int) -> Path:
    return self.trajectories_to_compare[i] if abs(i) < len(self) else None

  def __len__(self) -> int:
    return len(self.trajectories_to_compare)

  def get_num_steps(self) -> int:
    return len(self.base_at_steps)

  def calculate_metric(
                        self,
                        metric_function: Callable,
                      ) -> Dict:

    res = {}
    res["to_base"] = self.calculate_metric_to_base(metric_function=metric_function)
    res["to_base_future"] = self.calculate_metric_to_base_future(metric_function=metric_function)
    res["base_to_start"] = self.calculate_metric_to_first_base(metric_function=metric_function)
    res["traj_to_start"] = [self.calculate_metric_to_first(
                              traj_idx=i,
                              metric_function=metric_function
                            ) for i in range(len(self))]

    return res


  def calulate_all_metrics(
                        self,
                        metric_funcs: List[Tuple[Callable, str]],
                      ) -> Dict:

    return { name: self.calculate_metric(metric_function=func) for func, name in metric_funcs }

  def calculate_metric_to_base(self, metric_function: Callable) -> List[List[float]]:
    num_steps = self.get_num_steps()
    res = []
    for i in range(num_steps):
      mols_at_i = [mols[i] for mols in self.mols_to_compare_at_steps]
      metric_values = au.calculate_metric_one_to_many(
                                                      target=self.base_at_steps[i],
                                                      mols=mols_at_i,
                                                      align=True,
                                                      metric_function=metric_function
                                                    )
      res.append(metric_values)
  
    return res

  def calculate_metric_to_base_future(self, metric_function: Callable) -> List[List[float]]:
    num_steps = self.get_num_steps()
    res = []
    for i in range(num_steps):
      metric_values = [0 for _ in range(len(self))]
      for traj_idx, mols in enumerate(self.mols_to_compare_at_steps):
        metric_values[traj_idx] = min([au.calculate_metric_between_two_molecules(
                                              target=self.base_at_steps[step_idx],
                                              mol=mols[step_idx],
                                              align=True,
                                              metric_function=metric_function
                                            )
                                            for step_idx in range(i, num_steps)]
                                      )

      res.append(metric_values)
  
    return res

  def calculate_metric_to_first(
                                self,
                                traj_idx: int,
                                metric_function: Callable
                              ) -> List[float]:

    mols = self.mols_to_compare_at_steps[traj_idx]
    res = au.calculate_metric_one_to_many(
                                        target=mols[0],
                                        mols=mols,
                                        align=True,
                                        metric_function=metric_function
                                      )  
    return res

  def calculate_metric_to_first_base(
                                self,
                                metric_function: Callable
                              ) -> List[float]:
    mols = self.base_at_steps
    res = au.calculate_metric_one_to_many(
                                        target=mols[0],
                                        mols=mols,
                                        align=True,
                                        metric_function=metric_function
                                      )  
    return res


  def calculate_rmsd_to_base(self) -> List[List[float]]:
    num_steps = self.get_num_steps()
    rmsd_vals = []
    for i in range(num_steps):
      mols_at_i = [mols[i] for mols in self.mols_to_compare_at_steps]
      rmsds_to_target = au.calculate_rmsd_one_to_many(
                                                      target=self.base_at_steps[i],
                                                      mols=mols_at_i,
                                                      align=True
                                                    )
      rmsd_vals.append(rmsds_to_target)
  
    return rmsd_vals

  def calculate_rmsd_to_base_future(self) -> List[List[float]]:
    num_steps = self.get_num_steps()
    rmsd_vals = []
    for i in range(num_steps):
      rmsds_to_target = [0 for _ in range(len(self))]
      for traj_idx, mols in enumerate(self.mols_to_compare_at_steps):
        rmsds_to_target[traj_idx] = min([au.calculate_rmsd(
                                              target=self.base_at_steps[step_idx],
                                              mol=mols[step_idx],
                                              align=True
                                            )
                                            for step_idx in range(i, num_steps)]
                                      )

      rmsd_vals.append(rmsds_to_target)
  
    return rmsd_vals

  def calculate_rmsd_to_first(self, traj_idx: int) -> List[float]:
    mols = self.mols_to_compare_at_steps[traj_idx]
    res = au.calculate_rmsd_one_to_many(
                                        target=mols[0],
                                        mols=mols,
                                        align=True
                                      )  
    return res

  def calculate_rmsd_to_first_base(self) -> List[float]:
    mols = self.base_at_steps
    res = au.calculate_rmsd_one_to_many(
                                        target=mols[0],
                                        mols=mols,
                                        align=True
                                      )  
    return res
