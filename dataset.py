from collections import namedtuple
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Set, Tuple, Union


@dataclass
class Dataset():
  description: str
  file_paths: List[Path] = field(default_factory=list)
  names: List[str] = field(init=False)


  def __post_init__(self):
    unique_paths: List[Path] = []
    keys: Set[str] = set()
    for f in self.file_paths:
      if f.is_file():
        key = self.get_key(f)
        if key not in keys:
          keys.add(key)
          unique_paths.append(f)

    self.file_paths = unique_paths
    self.names = [f.name for f in unique_paths]


  def __getitem__(self, i: int) -> Path:
    return self.file_paths[i] if i >= 0 and i < len(self) else None


  def __len__(self) -> int:
    return len(self.file_paths)


  @staticmethod
  def create_with_names(description: str, paths_with_names: List[Tuple[Path, str]]):
    '''
      Give each filepath a unique name. That helps to track where the data is coming from.
    '''
    res = Dataset(description=description, file_paths=[])
    for file_path, name in paths_with_names:
      res.add(file_path=file_path, name=name)

    return res


  def add(self, file_path: Path, name: str=None) -> int:
    f_path = Path(file_path)
    if f_path.is_file() and self.is_unique(f_path):
      self.file_paths.append(f_path)
      if name and len(name.strip()) > 0:
        self.names.append(name)
      else:
        self.names.append(f_path.name)
      return len(self)

    return -1


  def add_many(self, file_paths: List[Path]):
    for x in file_paths:
      self.add(x)


  def clear(self):
    self.file_paths.clear()
    self.names.clear()


  def get_key(self, file_path: Path) -> str:
    return str(Path(file_path).absolute())


  def is_unique(self, file_path: Path) -> bool:
    key = self.get_key(file_path)
    for f in self.file_paths:
      if self.get_key(f) == key:
        return False
    return True

  # def merge_with(self, other_dataset: Dataset) -> Dataset:
  #   return merge_datasets(self, other_dataset)

  def remove(self, idx: int):
    if self[idx] != None:
      self.file_paths.pop(idx)
      self.names.pop(idx)


def merge_datasets(
                    dataset1: Dataset,
                    dataset2: Dataset,
                    description: str
                  ):

  res = Dataset.create_with_names(
                    description=description,
                    paths_with_names=[(p,n) for p,n in zip(dataset1.file_paths, dataset1.names)]
                  )

  for p,n in zip(dataset2.file_paths, dataset2.names):
    res.add(file_path=p, name=n)

  return res
