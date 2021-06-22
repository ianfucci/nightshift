from collections import namedtuple
import itertools
import logging
from operator import attrgetter
from typing import Dict, List, NamedTuple, Tuple

from nightshift import constants
from nightshift.bmrb import BMRBEntity

class Selector:
    '''Defines selector for any correlations in constants.SIMPLE_ATOMS'''
    def __init__(self, residues: List[str], atoms: Tuple[str], segment: Tuple[int] = None) -> None:
        self.residues = residues
        self.atoms = atoms
        self.segment = segment
    
    @property
    def selections(self) -> Dict[str, List[Tuple[str]]]:
        return {residue: [self.atoms] for residue in self.residues}

    def get_correlations(self, entity: BMRBEntity) -> List:
        residue_shifts = self._get_segment(entity)
        shift_table = []
        for seq_num, residue in residue_shifts.items():
            try:
                residue_type = residue[0].Comp_ID
                for correlation in self.selections[residue_type]:
                    # correlation has to be  first in the list comp to preserve order
                    selected_atoms = tuple(float(atom.Val)
                                           for spin in correlation
                                           for atom in residue 
                                           if atom.Atom_ID == spin)
                    if len(selected_atoms) == len(self.atoms):
                        shift_table.append([seq_num, residue_type, selected_atoms])
            except KeyError:
                # residue filtered out
                continue
        return shift_table    

    def _get_segment(self, entity: BMRBEntity) -> Dict[int, List[NamedTuple]]:
        # Group atoms by residue, each item is all assigned shifts for that residue
        residue_shifts = {int(seq_num) : list(group)
                         for seq_num, group 
                         in itertools.groupby(entity.shifts, key=attrgetter('Seq_ID'))
                         }

        # Apply segment filter
        sequence_min = min(residue_shifts.keys())
        sequence_max = max(residue_shifts.keys())
        start, stop = sequence_min, sequence_max
        if self.segment is not None:
            start, stop = self.segment
            if start < sequence_min:
                logging.warn(f'Start from segment parameter: {start} too low, defaulting to sequence start: residue {sequence_min}')
                start = sequence_min 
            if stop > sequence_max:
                logging.warn(f'Stop from segment parameter: {stop} too high, defaulting to sequence end: residue {sequence_max}')
                stop = sequence_max
            if start > sequence_max:
                logging.warn(f'Start from segment parameter: {stop} too high, ignoring')
                start = sequence_min
        return {seq_num : group
               for seq_num, group in residue_shifts.items()
               if start <= seq_num <= stop
               }

class AmideSelector(Selector):
    def __init__(self, residues: List[str], segment: Tuple[int] = None,
                *, sidechains: bool = False) -> None:
        super().__init__(residues, ('H', 'N'), segment)
        self.sidechains = sidechains

    @property
    def selections(self) -> Dict[str, List[Tuple[str]]]:
        residue_selections = super().selections
        if 'TRP' in self.residues:
            residue_selections['TRP'].append(('HE1', 'NE1'))
        if self.sidechains:
            if 'ASN' in self.residues:
                residue_selections['ASN'].append(('HD21','ND2'))
                residue_selections['ASN'].append(('HD22','ND2'))
            if 'GLN' in self.residues:
                residue_selections['GLN'].append(('HE21','NE2'))
                residue_selections['GLN'].append(('HE22','NE2'))
        return residue_selections

class MethylSelector(Selector):
    def __init__(self, residues: List[str], segment: Tuple[int] = None,
                *, proR: bool = False, proS: bool = False) -> None:
        methyl_residues = [residue 
                           for residue in residues
                           if residue in constants.METHYL_ATOMS.keys()]
        super().__init__(methyl_residues, ('HMETHYL', 'CMETHYL'), segment)
        self.proR = proR
        self.proS = proS
    
    @property
    def selections(self) -> Dict[str, List[Tuple[str]]]:
        if self.proR:
            methyl_selection = constants.METHYL_ATOMS_PROR
        elif self.proS:
            methyl_selection = constants.METHYL_ATOMS_PROS
        else:
            methyl_selection = constants.METHYL_ATOMS
        return {residue: atoms 
               for residue, atoms in methyl_selection.items()
               if residue in self.residues
               }

class AdvancedSelector(Selector):
    def __init__(self, residues: List[str], atoms: Tuple[str], segment: Tuple[int] = None,
                *, plus_minus: List[int], proR: bool = False, proS: bool = False) -> None:
        super().__init__(residues, atoms, segment)
        self.plus_minus = plus_minus
        self.proR = proR
        self.proS = proS

    @property
    def selections(self) -> Dict[str, List[Tuple[str]]]:
        residue_selections = {}
        for residue in self.residues:
            selected_atoms = []
            for atom in self.atoms:
                if atom in constants.SIMPLE_ATOMS:
                    selected_atoms.append([atom])
                elif residue in constants.METHYL_ATOMS.keys() and (atom == 'HMETHYL' or atom == 'CMETHYL'):
                    if self.proR:
                        methyl_selection = constants.METHYL_ATOMS_PROR
                    elif self.proS:
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
    
    def get_correlations(self, entity: BMRBEntity, *, label=None) -> List[NamedTuple]:
        # Intra-residue correlation same as the simple case
        if self.plus_minus == [0]*len(self.atoms):
            return super().get_correlations(entity)
        
        # Inter-residue correlation
        else:
            # Group atoms by residue, each item is all assigned shifts for that residue
            residue_shifts = self._get_segment(entity)
            shift_table = []
            # Filter out residues if sequence numbers for i+- cannot be found
            for ires in residue_shifts:
                window = tuple(ires + index for index in self.plus_minus)
                residues = []
                for res_num in window:
                    try:
                        residues.append(residue_shifts[res_num])
                    except KeyError:
                        continue
                
                # Get residue number, type and chemical shifts for each atom in group
                group = []
                for i, residue in enumerate(residues):
                    try:
                        # Have to filter it to prevent duplications in branched amino acid residues
                        correlations = [s[i] for s in self.selections[residue[0].Comp_ID]]
                        unique_correlations = sorted(set(correlations), key=correlations.index)

                        selected_atoms = []
                        for correlation in unique_correlations:
                            for atom in residue:
                                if atom.Atom_ID == correlation:
                                    selected_atoms.append((int(atom.Seq_ID), atom.Comp_ID, float(atom.Val)))
                        group.append(selected_atoms)
                    except KeyError:
                        # residue filtered out
                        continue

                # Reformat data and filter by those which have all atoms assigned
                for group_data in itertools.product(*group):
                    group_sequence_numbers, group_residue_types, group_shifts = zip(*group_data)
                    try:
                        if label is not None:
                            sequence_number = group_sequence_numbers[label]
                            residue_type = group_residue_types[label]
                        else:
                            ires_index = self.plus_minus.index(0)
                            sequence_number = group_sequence_numbers[ires_index]
                            residue_type = group_residue_types[ires_index]
                        if len(group_shifts) == len(self.atoms):
                            shift_table.append([sequence_number, residue_type, group_shifts])
                    except IndexError:
                        # not all atoms have assigned shifts
                        continue
        return shift_table