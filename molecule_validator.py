from dataclasses import dataclass, field
import numpy as np
from typing import Any, Dict, List, Tuple, Union

from ase import Atoms

import ase_utils as au
import features as ft
from features_constraints import Constraint, FeatureConstraint, MoleculeFeatureConstraints


def validate_structure(
                        mol: Atoms,
                        constraints: List[FeatureConstraint]
                      ) -> bool:

  return MoleculeFeatureConstraints(
              label="",
              molecule=mol,
              feature_constraints=constraints
            ).satisfied()



def target_mol_constraints() -> List[FeatureConstraint]:
  return [
    FeatureConstraint(
        label="",
        feature=ft.Angle("", 0, 0, 0),
        constraint=Constraint("", "<", 0)
      )


  ]

