from collections import namedtuple
import itertools
from operator import attrgetter
from typing import Dict, List, NamedTuple, Tuple

from nightshift import constants
from nightshift.bmrb import BMRBEntity

class Selector:
    '''Defines selector for any correlations in constants.SIMPLE_ATOMS'''
    def __init__(self, residues: List[str], atoms: Tuple[str]) -> None:
        self.residues = residues
        self.atoms = atoms
    
    def get_selections(self) -> Dict[str, List[Tuple[str]]]:
        return {residue: [self.atoms] for residue in self.residues}

    def get_correlations(self, entity: BMRBEntity, **kwargs) -> List[NamedTuple]:
        # Group atoms by residue with each item is all assigned shifts for that residue
        selections = self.get_selections(**kwargs) # allows subclasses to send kwargs
        residue_shifts = [list(g) for _, g in itertools.groupby(entity.shifts, key=attrgetter('Seq_ID'))]
        NightshiftRecord = namedtuple('NightshiftRecord', ['label'] + list(self.atoms))
        records = []
        for residue in residue_shifts:
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

class AmideSelector(Selector):
    def __init__(self, residues: List[str]) -> None:
        super().__init__(residues, ('H', 'N'))
    
    def get_selections(self, *, sidechains=False) -> Dict[str, List[Tuple[str]]]:
        # Potentially add a flag to show side chain GLN/ASN
        residue_selections = super().get_selections()
        if 'TRP' in self.residues:
            residue_selections['TRP'].append(('HE1', 'NE1'))
        if sidechains:
            if 'ASN' in self.residues:
                residue_selections['ASN'].append(('HD21','ND2'))
                residue_selections['ASN'].append(('HD22','ND2'))
            if 'GLN' in self.residues:
                residue_selections['GLN'].append(('HE21','NE2'))
                residue_selections['GLN'].append(('HE22','NE2'))
        return residue_selections

class MethylSelector(Selector):
    def __init__(self, residues: List[str]) -> None:
        methyl_residues = [residue 
                           for residue in residues
                           if residue in constants.METHYL_ATOMS.keys()]
        super().__init__(methyl_residues, ('HMETHYL', 'CMETHYL'))
    
    def get_selections(self, *, proR=False, proS=False) -> Dict[str, List[Tuple[str]]]:
        if proR:
            methyl_selection = constants.METHYL_ATOMS_PROR
        elif proS:
            methyl_selection = constants.METHYL_ATOMS_PROS
        else:
            methyl_selection = constants.METHYL_ATOMS
        return {residue: atoms 
                for residue, atoms in methyl_selection.items()
                if residue in self.residues}

class AdvancedSelector(Selector):
    def __init__(self, residues: List[str], atoms: Tuple[str], *, plus_minus: List[int]) -> None:
        super().__init__(residues, atoms)
        self.plus_minus = plus_minus

    def get_selections(self, *, proR=False, proS=False) -> Dict[str, List[Tuple[str]]]:
        residue_selections = {}
        for residue in self.residues:
            selected_atoms = []
            for atom in self.atoms:
                if atom in constants.SIMPLE_ATOMS:
                    selected_atoms.append([atom])
                elif residue in constants.METHYL_ATOMS.keys() and (atom == 'HMETHYL' or atom == 'CMETHYL'):
                    if proR:
                        methyl_selection = constants.METHYL_ATOMS_PROR
                    elif proS:
                        methyl_selection = constants.METHYL_ATOMS_PROS
                    else:
                        methyl_selection = constants.METHYL_ATOMS
                    # first atom is always proton, second is carbon in METHYL_ATOMS dicts
                    index = 0 if atom == 'HMETHYL' else 1
                    selected_atoms.append([methyl_atoms[index] 
                                          for methyl_atoms in methyl_selection[residue]])
                else:
                    selected_atoms.append([residue_atom 
                                           for residue_atom in constants.RESIDUE_ATOMS[residue] 
                                           if residue_atom.startswith(atom)])
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