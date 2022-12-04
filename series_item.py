from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union


@dataclass
class SeriesItem():
  name: str
  labels: List[str] = field(default_factory=list)
  values: List[float] = field(default_factory=list)
