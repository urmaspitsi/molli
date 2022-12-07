from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple, Union
import numpy as np

from ase import Atoms
import features as ft
import metrics as ms


@dataclass
class Constraint():
  label: str
  condition: str
  value: float
  abs_tolerance: float = 0.0000001

  def satisfied(self, test_with_value: float) -> bool:
    if self.condition == "<":
      return test_with_value < self.value
    elif self.condition == ">":
      return test_with_value > self.value
    elif self.condition == "<=":
      return test_with_value <= self.value
    elif self.condition == ">=":
      return test_with_value >= self.value
    elif self.condition == "=":
      return abs(self.value - test_with_value) < self.abs_tolerance


@dataclass
class FeatureConstraint():
  label: str
  feature: ft.Feature
  constraint: Constraint


@dataclass
class MoleculeFeatureConstraints():
  label: str
  molecule: Atoms
  feature_constraints: List[FeatureConstraint]

  def __call__(self, *args: Any, **kwds: Any) -> bool:
    return self.satisfied()

  def satisfied(self) -> bool:
    mol = self.molecule

    for fc in self.feature_constraints:
      if not fc.constraint.satisfied(test_with_value=fc.feature.calculate_value(mol)):
        return False

    return True
