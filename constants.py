from ase import units as ase_units


atomic_numbers_to_elements = {
    1:'H',2:'He',3:'Li',4:'Be',5:'B',6:'C',7:'N',8:'O',9:'F',10:'Ne',\
    11:'Na',12:'Mg',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',18:'Ar',19:'K',\
    20:'Ca',21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',27:'Co',28:'Ni',\
    29:'Cu',30:'Zn',31:'Ga',32:'Ge',33:'As',34:'Se',35:'Br',36:'Kr',37:'Rb',\
    38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',44:'Ru',45:'Rh',46:'Pd',\
    47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',53:'I',54:'Xe',55:'Cs',\
    56:'Ba',57:'La',58:'Ce',59:'Pr',60:'Nd',61:'Pm',62:'Sm',63:'Eu',64:'Gd',\
    65:'Tb',66:'Dy',67:'Ho',68:'Er',69:'Tm',70:'Yb',71:'Lu',72:'Hf',73:'Ta',\
    74:'W',75:'Re',76:'Os',77:'Ir',78:'Pt',79:'Au',80:'Hg',81:'Tl',82:'Pb',\
    83:'Bi',84:'Po',85:'At',86:'Rn',87:'Fr',88:'Ra',89:'Ac',90:'Th',91:'Pa',\
    92:'U',93:'Np',94:'Pu',95:'Am',96:'Cm',97:'Bk',98:'Cf',99:'Es',100:'Fm',
  }

hartree_in_kcal_per_mol = ase_units.Hartree * ase_units.mol / ase_units.kcal

hartree_in_kJ_per_mol = ase_units.Hartree * ase_units.mol / ase_units.kJ


def convert_elements_to_numbers():
  elem2nr = {}
  for k,v in atomic_numbers_to_elements.items():
    elem2nr[v.lower()] = k
  return elem2nr

elements_to_atomic_numbers = convert_elements_to_numbers()


def get_atomic_number(element: str) -> int:
  return elements_to_atomic_numbers[element.lower()]

def get_element(atomic_nr: int) -> str:
  return atomic_numbers_to_elements[atomic_nr]
