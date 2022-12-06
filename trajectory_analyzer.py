from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, NamedTuple, Set, Tuple, Union
import numpy as np

from ase import Atoms

import ase_utils as au
from file_source import MultiItemFileSource
from series_item import SeriesItem

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
    vs_base = self.calculate_metric_to_base(metric_function=metric_function)
    vs_base_future = self.calculate_metric_to_base_future(metric_function=metric_function)
    vs_base_last_step = self.calculate_metric_to_base_last_step(metric_function=metric_function)

    res["base_name"] = self.base_trajectory.name
    res["names_to_compare"] = [x.name for x in self.trajectories_to_compare]
    res["to_base"] = np.array(vs_base).transpose()
    res["to_base_future"] = np.array(vs_base_future).transpose()
    res["to_base_last_step"] = np.array(vs_base_last_step).transpose()
    res["base_to_start"] = self.calculate_metric_to_first_base(metric_function=metric_function)
    res["base_to_last"] = self.calculate_metric_to_last_base(metric_function=metric_function)
    res["traj_to_start"] = [self.calculate_metric_to_first(
                              traj_idx=i,
                              metric_function=metric_function
                            ) for i in range(len(self))]
    res["base_incremental"] = self.calculate_metric_incremental_base(metric_function=metric_function)
    res["traj_incremental"] = [self.calculate_metric_incremental(
                              traj_idx=i,
                              metric_function=metric_function
                            ) for i in range(len(self))]

    return res


  def calulate_all_metrics(
                        self,
                        metric_funcs: List[Tuple[Callable, str]],
                      ) -> Dict:

    return { name: self.calculate_metric(metric_function=func) for func, name in metric_funcs }


  def calculate_metric_incremental(
                                self,
                                traj_idx: int,
                                metric_function: Callable
                              ) -> List[float]:

    mols = self.mols_to_compare_at_steps[traj_idx]
    res = au.calculate_metric_many_to_many(
                                        targets=mols[:-1],
                                        mols=mols[1:],
                                        align=True,
                                        metric_function=metric_function
                                      )  
    return res


  def calculate_metric_incremental_base(
                                self,
                                metric_function: Callable
                              ) -> List[float]:

    mols = self.base_at_steps
    res = au.calculate_metric_many_to_many(
                                        targets=mols[:-1],
                                        mols=mols[1:],
                                        align=True,
                                        metric_function=metric_function
                                      )  
    return res


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


  def calculate_metric_to_base_last_step(self, metric_function: Callable) -> List[List[float]]:
    base_at_last_step = self.base_at_steps[-1]
    num_steps = self.get_num_steps()
    res = []
    for i in range(num_steps):
      mols_at_i = [mols[i] for mols in self.mols_to_compare_at_steps]
      metric_values = au.calculate_metric_one_to_many(
                                                      target=base_at_last_step,
                                                      mols=mols_at_i,
                                                      align=True,
                                                      metric_function=metric_function
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


  def calculate_metric_to_last(
                                self,
                                traj_idx: int,
                                metric_function: Callable
                              ) -> List[float]:

    mols = self.mols_to_compare_at_steps[traj_idx]
    res = au.calculate_metric_one_to_many(
                                        target=mols[-1],
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


  def calculate_metric_to_last_base(
                                self,
                                metric_function: Callable
                              ) -> List[float]:
    mols = self.base_at_steps
    res = au.calculate_metric_one_to_many(
                                        target=mols[-1],
                                        mols=mols,
                                        align=True,
                                        metric_function=metric_function
                                      )  
    return res
