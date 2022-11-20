import features as FT


mol24_forbidden_bonds = [
  FT.Distance(label="7-2", atom_idx1=7, atom_idx2=2),
  FT.Distance(label="7-12", atom_idx1=7, atom_idx2=12),
  FT.Distance(label="61-45", atom_idx1=61, atom_idx2=45),
  FT.Distance(label="61-65", atom_idx1=61, atom_idx2=65),
]

mol24_constraints6_angles = [
  FT.Angle(label="7-6-4", atom_idx1=7, atom_idx2=6, atom_idx3=4),
  FT.Angle("7-20-18", 7, 20, 18),
  FT.Angle("7-31-29", 7, 31, 29),
  FT.Angle("7-60-53", 7, 60, 53),
  FT.Angle("61-56-48", 61, 56, 48),
  FT.Angle("61-59-54", 61, 59, 54),
  FT.Angle("61-60-53", 61, 60, 53),
  FT.Angle("61-68-72", 61, 68, 72),
  FT.Angle("7-8-10", 7, 8, 10),
  FT.Angle("7-8-61", 7, 8, 61),
  FT.Angle("7-60-61", 7, 60, 61),
]

