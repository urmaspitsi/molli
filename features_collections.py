import features as FT


mol24_forbidden_bonds = [
  FT.Distance(label="distance Ti(7)-O(2)", atom_idx1=7-1, atom_idx2=2-1),
  FT.Distance(label="distance Ti(7)-O(12)", atom_idx1=7-1, atom_idx2=12-1),
  FT.Distance(label="distance Ti(61)-O(45)", atom_idx1=61-1, atom_idx2=45-1),
  FT.Distance(label="distance Ti(61)-O(65)", atom_idx1=61-1, atom_idx2=65-1),
]

mol24_constraints6_angles = [
  FT.Angle(label="angle Ti(7)-O(6)-C(4)", atom_idx1=7-1, atom_idx2=6-1, atom_idx3=4-1),
  FT.Angle("angle Ti(7)-O(20)-C(18)", 7-1, 20-1, 18-1),
  FT.Angle("angle Ti(7)-O(31)-C(29)", 7-1, 31-1, 29-1),
  FT.Angle("angle Ti(7)-O(60)-C(53)", 7-1, 60-1, 53-1),
  FT.Angle("angle Ti(61)-O(56)-C(48)", 61-1, 56-1, 48-1),
  FT.Angle("angle Ti(61)-O(59)-C(54)", 61-1, 59-1, 54-1),
  FT.Angle("angle Ti(61)-O(60)-C(53)", 61-1, 60-1, 53-1),
  FT.Angle("angle Ti(61)-O(68)-C(72)", 61-1, 68-1, 72-1),
  FT.Angle("angle Ti(7)-O(8)-C(10)", 7-1, 8-1, 10-1),
  FT.Angle("angle Ti(7)-O(8)-Ti(61)", 7-1, 8-1, 61-1),
  FT.Angle("angle Ti(7)-O(60)-Ti(61)", 7-1, 60-1, 61-1),
]

