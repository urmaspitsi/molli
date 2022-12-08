from dataclasses import dataclass, field
from typing import Any, Dict, List, Tuple, Union

# @dataclass
# class Bonds():


@dataclass
class TargetMoleculeBonds():
  '''
    created with ase.
    source: source: https://wiki.fysik.dtu.dk/ase/ase/geometry.html
    from ase.geometry.analysis import Analysis
    mol = ase_io.read("..xyz")
    ana = Analysis(mol)
    ana.get_bonds("C", "C", unique=True)
  '''


  CO = [
      (2, 0),
      (2, 1),
      (3, 5),
      (9, 7),
      (10, 11),
      (10, 12),
      (17, 19),
      (28, 30),
      (36, 37),
      (41, 37),
      (41, 44),
      (47, 55),
      (52, 59),
      (53, 58),
      (61, 64),
      (61, 65),
      (62, 65),
      (71, 67),
      (78, 12),
      (82, 0)
      ]

  CC = [
    (2, 3),
    (3, 9),
    (9, 10),
    (15, 17),
    (17, 21),
    (26, 28),
    (28, 32),
    (41, 52),
    (45, 47),
    (47, 50),
    (52, 53),
    (53, 61),
    (70, 71),
    (71, 73)
    ]

  CH = [
    (3, 4),
    (9, 8),
    (15, 13),
    (15, 14),
    (15, 18),
    (17, 16),
    (21, 20),
    (21, 22),
    (21, 23),
    (26, 24),
    (26, 25),
    (26, 29),
    (28, 27),
    (32, 31),
    (32, 33),
    (32, 34),
    (36, 35),
    (36, 38),
    (36, 39),
    (45, 40),
    (45, 42),
    (45, 49),
    (47, 43),
    (50, 46),
    (50, 54),
    (50, 57),
    (52, 51),
    (53, 48),
    (62, 56),
    (62, 63),
    (62, 66),
    (70, 68),
    (70, 69),
    (70, 74),
    (71, 75),
    (73, 72),
    (73, 76),
    (73, 77),
    (78, 79),
    (78, 80),
    (78, 81),
    (82, 83),
    (82, 84),
    (82, 85)
    ]

  bonds_distances = {
    (2,3): 1.547103657665833,
    (3,9): 1.5894456376419421,
    (9,10): 1.546822288655035,
    (15,17): 1.5356032980568255,
    (17,21): 1.5347359309053787,
    (26,28): 1.534118534196755,
    (28,32): 1.5337650488419012,
    (41,52): 1.5355396336900589,
    (45,47): 1.5340929182086722,
    (47,50): 1.5345500413639173,
    (52,53): 1.5745448621423905,
    (53,61): 1.5442633948977098,
    (70,71): 1.5344868420400353,
    (71,73): 1.533219391773402,
    (3,4): 1.1138821499660545,
    (8,9): 1.1147429993541111,
    (13,15): 1.1109931689632482,
    (14,15): 1.1117729366246507,
    (15,18): 1.1129948796333253,
    (16,17): 1.1189068301757745,
    (20,21): 1.1115617320009719,
    (21,22): 1.1116485572342547,
    (21,23): 1.112768121946347,
    (24,26): 1.1118906755257012,
    (25,26): 1.1125975872210043,
    (26,29): 1.1119184193042222,
    (27,28): 1.1169407097281392,
    (31,32): 1.1128204211915775,
    (32,33): 1.1100831087648346,
    (32,34): 1.112588553633822,
    (35,36): 1.1064985486271548,
    (36,38): 1.1104563954577413,
    (36,39): 1.1101082863221048,
    (40,45): 1.112937967554347,
    (42,45): 1.1110876799262963,
    (45,49): 1.1132531567056076,
    (43,47): 1.1157168902091605,
    (46,50): 1.112996050577,
    (50,54): 1.1116985103457684,
    (50,57): 1.1129139058215596,
    (51,52): 1.1092533610956519,
    (48,53): 1.121130039515934,
    (56,62): 1.1101532095715436,
    (62,63): 1.1066072678281131,
    (62,66): 1.1103661811226961,
    (68,70): 1.1123679582408872,
    (69,70): 1.1125454969595625,
    (70,74): 1.1127837529259672,
    (71,75): 1.1174792497299448,
    (72,73): 1.113335061906343,
    (73,76): 1.1106301512047112,
    (73,77): 1.113243797684047,
    (78,79): 1.1105758710614957,
    (78,80): 1.1118084106036434,
    (78,81): 1.1077124097025361,
    (82,83): 1.106371656510144,
    (82,84): 1.109770054329274,
    (82,85): 1.110147417556335,
    (0,2): 1.3345357174159858,
    (1,2): 1.231020339583388,
    (3,5): 1.3850773860777597,
    (7,9): 1.3920563433514466,
    (10,11): 1.2133539583518074,
    (10,12): 1.3617756741159681,
    (17,19): 1.4118194556635066,
    (28,30): 1.4164519906812938,
    (36,37): 1.4451431472369787,
    (37,41): 1.3361156056958547,
    (41,44): 1.2304660062151251,
    (47,55): 1.4174629172803783,
    (52,59): 1.401873379823941,
    (53,58): 1.381021849028465,
    (61,64): 1.216136370083964,
    (61,65): 1.3529672579800296,
    (62,65): 1.4424313938655802,
    (67,71): 1.416957859701198,
    (12,78): 1.435528289307111,
    (0,82): 1.442119571821283
}