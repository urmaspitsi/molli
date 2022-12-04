from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Set, Tuple, Union


@dataclass
class FileSource():
  file_path: Path
  name: str


@dataclass
class MultiItemFileSource():
  file_path: Path
  name: str
  item_idxs: List[int] = field(default_factory=list)

