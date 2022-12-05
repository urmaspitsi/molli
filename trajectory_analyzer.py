from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Set, Tuple, Union
from ase import Atoms

import ase_utils as au
from file_source import MultiItemFileSource

@dataclass
class TrajectoryAnalyzer():
  description: str
  base_trajectory: MultiItemFileSource
  trajectories_to_compare: List[MultiItemFileSource] = field(default_factory=list)
  base_at_steps: List[Atoms] = field(init=False)
  mols_to_compare_at_steps: List[List[Atoms]] = field(init=False)

  def __post_init__(self):
    self.base_at_steps = au.extract_mols_from_xyz_files(sources=[self.base_trajectory])[0]
    self.mols_to_compare_at_steps = au.extract_mols_from_xyz_files(sources=self.trajectories_to_compare)

  def __getitem__(self, i: int) -> Path:
    return self.trajectories_to_compare[i] if abs(i) < len(self) else None

  def __len__(self) -> int:
    return len(self.trajectories_to_compare)

  def calculate_rmsd_to_base(self) -> List[List[float]]:
    num_step_idxs = len(self.base_trajectory.item_idxs)
    rmsd_vals = []
    for i in range(num_step_idxs):
      mols_at_i = [mols[i] for mols in self.mols_to_compare_at_steps]
      rmsds_to_target = au.calculate_rmsd_one_to_many(
                                                      target=self.base_at_steps[i],
                                                      mols=mols_at_i,
                                                      align=True
                                                    )
      rmsd_vals.append(rmsds_to_target)
  
    return rmsd_vals

  def calculate_rmsd_to_base_future(self) -> List[List[float]]:
    num_step_idxs = len(self.base_trajectory.item_idxs)
    rmsd_vals = []
    for i in range(num_step_idxs):
      rmsds_to_target = [0 for x in range(len(self))]
      for traj_idx, mols in enumerate(self.mols_to_compare_at_steps):
        rmsds_to_target[traj_idx] = min([au.calculate_rmsd(
                                              target=self.base_at_steps[step_idx],
                                              mol=mols[step_idx],
                                              align=True
                                            )
                                            for step_idx in range(i, num_step_idxs)]
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