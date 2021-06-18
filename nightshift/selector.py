from collections import namedtuple
import itertools
from operator import attrgetter
from typing import Dict, List, NamedTuple, Tuple

from nightshift import constants
from nightshift.bmrb import BMRBEntity

class Selector:
    '''Custom correlations'''
    def __init__(self, residues: List[str], atoms: Tuple[str], *, plus_minus: List[int]) -> None:
        self.residues = residues
        self.atoms = atoms
        self.plus_minus = plus_minus
    
    def get_selections(self, *, proR=False, proS=False) -> Dict[str, List[Tuple[str]]]:
        residue_selections = {}
        for residue in self.residues:
            selected_atoms = []
            for atom in self.atoms:
                if atom in constants.SIMPLE_ATOMS:
                    selected_atoms.append([atom])
                elif residue in constants.METHYL_ATOMS.keys() and atom == 'HMETHYL':
                    pairs = constants.METHYL_ATOMS[residue]
                    if residue in {'LEU', 'VAL'}:
                        if proR:
                            selected_atoms.append([p[0] for p in pairs if not p[1].endswith('2')])
                        elif proS:
                            selected_atoms.append([p[0] for p in pairs if not p[1].endswith('1')])
                        else:
                            # Not proR or proS filtered
                            selected_atoms.append([p[0] for p in pairs])
                    else:
                        # Not Leu or Val
                        selected_atoms.append([p[0] for p in pairs])
                elif residue in constants.METHYL_ATOMS.keys() and atom == 'CMETHYL':
                    pairs = constants.METHYL_ATOMS[residue]
                    if residue in {'LEU', 'VAL'}:
                        if proR:
                            selected_atoms.append([p[1] for p in pairs if not p[1].endswith('2')])
                        elif proS:
                            selected_atoms.append([p[1] for p in pairs if not p[1].endswith('1')])
                        else:
                            # Not proR or proS filtered
                            selected_atoms.append([p[1] for p in pairs])
                    else:
                        # Not Leu or Val
                        selected_atoms.append([p[1] for p in pairs])
                else:
                    selected_atoms.append([res_atom 
                                           for res_atom in constants.RESIDUE_ATOMS[residue] 
                                           if res_atom.startswith(atom)])
            residue_selections[residue] = list(itertools.product(*selected_atoms))
        return residue_selections   

    def get_correlations(self, entity: BMRBEntity, *, label=None, **kwargs) -> List[NamedTuple]:
        # Group atoms by residue with each item is all assigned shifts for that residue
        selections = self.get_selections(**kwargs) # have to account for proR proS
        residue_shifts = {int(k) : list(g) for k, g in itertools.groupby(entity.shifts, key=attrgetter('Seq_ID'))}
        NightshiftRecord = namedtuple('NightshiftRecord', ['label'] + list(self.atoms))
        records = []
        if self.plus_minus == [0]*len(self.atoms):
            # Intra-residue
            for residue in residue_shifts.values():
                try:
                    residue_type = residue[0].Comp_ID
                    label = residue[0].Comp_ID + residue[0].Seq_ID
                    correlations = selections[residue_type]
                    for correlation in correlations:
                        # correlation has to be  first in the list comp to preserve order
                        selected_atoms = [float(atom.Val)
                                          for spin in correlation
                                          for atom in residue 
                                          if atom.Atom_ID == spin]
                        records.append(NightshiftRecord(label, *selected_atoms))
                except (KeyError, TypeError):
                    # KeyError: residue filtered out
                    # TypeError: not all atoms in correlation have assigned shifts
                    continue
            return records
              
        else:
            # Inter-residue correlation
            for ires in residue_shifts:
                window = tuple(ires + index for index in self.plus_minus)
                residues = []
                for res_num in window:
                    try:
                        residues.append(residue_shifts[res_num])
                    except KeyError:
                        # Sequence with i+- plus_minus not found
                        continue
                group = []
                for i, residue in enumerate(residues):
                    try:
                        residue_type = residue[0].Comp_ID

                        # For branched amino acids multiple tuples are in each list
                        # Have to filter it to prevent duplications in strange places.
                        correlations = [s[i] for s in selections[residue_type] if s]
                        unique_correlations = sorted(set(correlations), key=correlations.index)

                        selected_atoms = []
                        for correlation in unique_correlations:
                            for atom in residue:
                                if atom.Atom_ID == correlation:
                                    selected_atoms.append((atom.Comp_ID + atom.Seq_ID, float(atom.Val)))
                        group.append(selected_atoms)
                    except KeyError:
                        # KeyError: residue filtered out
                        continue
                for item in itertools.product(*group):
                    residue_labels = [field[0] for field in item]
                    residue_values = [field[1] for field in item]
                    try:
                        if label is not None:
                            record_label = residue_labels[label]
                        else:
                            ires_index = self.plus_minus.index(0)
                            record_label = residue_labels[ires_index]
                        records.append(NightshiftRecord(record_label, *residue_values))
                    except (IndexError, TypeError):
                        # not all atoms have assigned shifts
                        continue
        return records
                    

class AmideSelector(Selector):
    def __init__(self, residues: List[str]) -> None:
        super().__init__(residues, ('H', 'N'), plus_minus=[0,0])
    
    def get_selections(self) -> Dict[str, List[Tuple[str]]]:
        # Potentially add a flag to show side chain GLN/ASN
        residue_selections = super().get_selections()
        if 'TRP' in self.residues:
            residue_selections['TRP'].append(('HE1', 'NE1'))
        return residue_selections
    
    def get_correlations(self, entity: BMRBEntity, *, label, **kwargs) -> List[NamedTuple]:
        return super().get_correlations(entity, label=label)

class MethylSelector(Selector):
    def __init__(self, residues: List[str]) -> None:
        methyl_residues = [residue 
                           for residue in residues
                           if residue in constants.METHYL_ATOMS.keys()]
        super().__init__(methyl_residues, ('HMETHYL', 'CMETHYL'), plus_minus=[0,0])
    
    def get_selections(self, *, proR=False, proS=False) -> Dict[str, List[Tuple[str]]]:
        # Need to do proR, proS
        methyl_selection = constants.METHYL_ATOMS
        # Filter out prochiral atoms when given proR/proS flags
        for residue, pairs in methyl_selection.items():
            if residue in {'LEU', 'VAL'}:
                if proR:
                    methyl_selection[residue] = [p for p in pairs if not p[1].endswith('2')]
                elif proS:
                    methyl_selection[residue] = [p for p in pairs if not p[1].endswith('1')]
        return {residue: atoms 
                for residue, atoms in methyl_selection.items()
                if residue in self.residues}