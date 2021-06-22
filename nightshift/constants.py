from typing import Dict, List, Set, Tuple

ONE_LETTER_TO_THREE_LETTER: Dict[str, str] = {
    'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU',
    'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'M':'MET',
    'K':'LYS', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP',
    'Y':'TYR', 'V':'VAL',
    }

SIMPLE_ATOMS: Set[str] = {'H', 'HA', 'N', 'C', 'CA', 'CB'} # atoms without numbers afterward

AMINO_ACID_WEIGHTS: Dict[str, float]= {
    'A': 89.0932, 'C': 121.1582, 'D': 133.1027, 'E': 147.1293, 'F': 165.1891, 
    'G': 75.0666, 'H': 155.1546, 'I': 131.1729, 'K': 146.1876, 'L': 131.1729,
    'M': 149.211, 'N': 132.1179, 'O': 255.3134, 'P': 115.1305, 'Q': 146.1445,
    'R': 174.201, 'S': 105.0926, 'T': 119.1192, 'U': 168.0532, 'V': 117.1463,
    'W': 204.225, 'Y': 181.1885,
    }

RESIDUE_ATOMS: Dict[str, Tuple[str]]= {
    'ALA': ('H', 'HA', 'HB1', 'HB2', 'HB3', 'C', 'CA', 'CB', 'N'), 
    'ARG': ('H', 'HA', 'HB2', 'HB3', 'HD2', 'HD3', 'HE', 'HG2', 'HG3', 'HH11', 'HH12', 'HH21', 'HH22', 'C', 'CA', 'CB', 'CD', 'CG', 'CZ', 'N', 'NE', 'NH1', 'NH2'), 
    'ASN': ('H', 'HA', 'HB2', 'HB3', 'HD21', 'HD22', 'C', 'CA', 'CB', 'CG', 'N', 'ND2'), 
    'ASP': ('H', 'HA', 'HB2', 'HB3', 'HD2', 'C', 'CA', 'CB', 'CG', 'N'), 
    'CYS': ('H', 'HA', 'HB2', 'HB3', 'HG', 'C', 'CA', 'CB', 'N'), 
    'GLN': ('H', 'HA', 'HB2', 'HB3', 'HE21', 'HE22', 'HG2', 'HG3', 'C', 'CA', 'CB', 'CD', 'CG', 'N', 'NE2'), 
    'GLU': ('H', 'HA', 'HB2', 'HB3', 'HE2', 'HG2', 'HG3', 'C', 'CA', 'CB', 'CD', 'CG', 'N'), 
    'GLY': ('H', 'H1', 'HA2', 'HA3', 'C', 'CA', 'N'), 
    'HIS': ('H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2', 'C', 'CA', 'CB', 'CD2', 'CE1', 'CG', 'N', 'ND1', 'NE2'), 
    'ILE': ('H', 'HA', 'HB', 'HG12', 'HG13', 'HD1', 'HD2', 'HD3', 'HG21', 'HG22', 'HG23', 'C', 'CA', 'CB', 'CD1', 'CG1', 'CG2', 'N'), 
    'LEU': ('H', 'HA', 'HB2', 'HB3', 'HG', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23', 'C', 'CA', 'CB', 'CD1', 'CD2', 'CG', 'N'), 
    'LYS': ('H', 'HA', 'HB2', 'HB3', 'HD2', 'HD3', 'HE2', 'HE3', 'HG2', 'HG3', 'C', 'CA', 'CB', 'CD', 'CE', 'CG', 'N', 'NZ', 'QZ'), 
    'MET': ('H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HE1', 'HE2', 'HE3', 'C', 'CA', 'CB', 'CE', 'CG', 'N'), 
    'PHE': ('H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ', 'C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'N'), 
    'PRO': ('H', 'H2', 'HA', 'HB2', 'HB3', 'HD2', 'HD3', 'HG2', 'HG3', 'C', 'CA', 'CB', 'CD', 'CG', 'N'), 
    'SER': ('H', 'HA', 'HB2', 'HB3', 'HG', 'C', 'CA', 'CB', 'N'), 
    'THR': ('H', 'HA', 'HB', 'HG1', 'HG21', 'HG22', 'HG23', 'C', 'CA', 'CB', 'CG2', 'N'), 
    'TRP': ('H', 'HA', 'HB2', 'HB3', 'HD1', 'HE1', 'HE3', 'HH2', 'HZ2', 'HZ3', 'C', 'CA', 'CB', 'CD1', 'CD2', 'CE2', 'CE3', 'CG', 'CH2', 'CZ2', 'CZ3', 'N', 'NE1'), 
    'TYR': ('H', 'HA', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1', 'HE2', 'HH', 'C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'N'), 
    'VAL': ('H', 'HA', 'HB', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23', 'C', 'CA', 'CB', 'CG1', 'CG2', 'N'),
}

METHYL_ATOMS: Dict[str, List[Tuple[str]]] = {'ILE':[('HD11','CD1')],
                'LEU':[('HD11','CD1'), ('HD21','CD2')],
                'VAL':[('HG11','CG1'), ('HG21','CG2')],
                'MET':[('HE1','CE')],
                'ALA':[('HB1','CB')],
                'THR':[('HG21','CG2')],
                }

METHYL_ATOMS_PROR: Dict[str, List[Tuple[str]]] = {'ILE':[('HD11','CD1')],
                'LEU':[('HD11','CD1')],
                'VAL':[('HG11','CG1')],
                'MET':[('HE1','CE')],
                'ALA':[('HB1','CB')],
                'THR':[('HG21','CG2')],
                }

METHYL_ATOMS_PROS: Dict[str, List[Tuple[str]]] = {'ILE':[('HD11','CD1')],
                'LEU':[('HD21','CD2')],
                'VAL':[('HG21','CG2')],
                'MET':[('HE1','CE')],
                'ALA':[('HB1','CB')],
                'THR':[('HG21','CG2')],
                }