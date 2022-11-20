import features as FT


mol24_forbidden_bonds = [
  FT.Distance(label="7-2", atom_idx1=7-1, atom_idx2=2-1),
  FT.Distance(label="7-12", atom_idx1=7-1, atom_idx2=12-1),
  FT.Distance(label="61-45", atom_idx1=61-1, atom_idx2=45-1),
  FT.Distance(label="61-65", atom_idx1=61-1, atom_idx2=65-1),
]

mol24_constraints6_angles = [
  FT.Angle(label="7-6-4", atom_idx1=7-1, atom_idx2=6-1, atom_idx3=4-1),
  FT.Angle("7-20-18", 7-1, 20-1, 18-1),
  FT.Angle("7-31-29", 7-1, 31-1, 29-1),
  FT.Angle("7-60-53", 7-1, 60-1, 53-1),
  FT.Angle("61-56-48", 61-1, 56-1, 48-1),
  FT.Angle("61-59-54", 61-1, 59-1, 54-1),
  FT.Angle("61-60-53", 61-1, 60-1, 53-1),
  FT.Angle("61-68-72", 61-1, 68-1, 72-1),
  FT.Angle("7-8-10", 7-1, 8-1, 10-1),
  FT.Angle("7-8-61", 7-1, 8-1, 61-1),
  FT.Angle("7-60-61", 7-1, 60-1, 61-1),
]

