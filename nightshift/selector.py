import itertools
import logging
from operator import attrgetter
import re
from typing import Dict, List, NamedTuple, Tuple

from nightshift import constants

Selections =  Dict[str, List[Tuple[str]]]
Correlations = List[Tuple[int,str,Tuple[float]]]
ResidueShifts = Dict[int, List[NamedTuple]]

class Selector:
    '''Defines selector for any correlations in constants.SIMPLE_ATOMS'''
    def __init__(self, 
                residues: List[str] = None,
                atoms: Tuple[str] = None,
                segment: Tuple[int] = None) \
                -> None:

        self.residues = self._filter_residues(residues)
        self.atoms = atoms
        self.segment = segment
    
    @property
    def selections(self) -> Selections:
        return {residue: [self.atoms] for residue in self.residues}

    def get_correlations(self, entity_shifts: List[NamedTuple]) -> Correlations:
        residue_shifts = self._get_residue_shifts(entity_shifts)
        shift_table = []
        for seq_num, residue in residue_shifts.items():
            residue_type = residue[0].Comp_ID
            selections = self.selections.get(residue_type)
            if selections is not None:
                for correlation in selections:
                    selected_atoms = []
                    # correlation has to be outer loop to preserve order
                    for spin in correlation:
                        for atom in residue:
                            if atom.Atom_ID == spin:
                                selected_atoms.append(float(atom.Val))
                    if len(selected_atoms) == len(self.atoms):
                        shift_table.append((seq_num, residue_type, tuple(selected_atoms)))
        return shift_table
    
    @staticmethod
    def _filter_residues(residues: List[str]) -> List[str]:
        if residues is not None:
            # Warn for incorrect 1-letter codes
            bad_codes = set(residues.upper()) - constants.ONE_LETTER_TO_THREE_LETTER.keys()
            if bad_codes:
                bad_code_string = ','.join(bad_codes)
                logging.warn(f'{bad_code_string} not valid 1-letter code(s)')
            # Remove bad codes and ignore
            residues = [constants.ONE_LETTER_TO_THREE_LETTER.get(residue.upper()) 
                        for residue in residues if residue not in bad_codes]
        else:
            # If no residue filter is specified use all residues
            residues = list(constants.ONE_LETTER_TO_THREE_LETTER.values())
        return residues

    def _get_residue_shifts(self, entity_shifts: List[NamedTuple]) -> ResidueShifts:
        # Group atoms by residue, each item is all assigned shifts for that residue
        residue_shifts = {int(seq_num) : list(group)
                         for seq_num, group 
                         in itertools.groupby(entity_shifts, key=attrgetter('Seq_ID'))
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
    def __init__(self,
                residues: List[str] = None, 
                segment: Tuple[int] = None,
                *,
                sidechains: bool = False) \
                -> None:

        super().__init__(residues, ('H', 'N'), segment)
        self.sidechains = sidechains

    @property
    def selections(self) -> Selections:
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
    def __init__(self,
                residues: List[str] = None,
                segment: Tuple[int] = None,
                *,
                proR: bool = False,
                proS: bool = False) \
                -> None:

        # Warn if non MILVAT residues are used with the --methyl flag
        super().__init__(residues, ('HMETHYL', 'CMETHYL'), segment)
        if residues is not None:
            non_milvat = set(self.residues) - constants.METHYL_ATOMS.keys()
            if non_milvat:
                res_string = ','.join(sorted(non_milvat, key=self.residues.index))
                logging.warn(f'residues other than MILVAT: ({res_string}) are ignored')

        self.residues = [residue 
                        for residue in constants.METHYL_ATOMS.keys()
                        if residue in self.residues]
        self.proR = proR
        self.proS = proS
    
    @property
    def selections(self) -> Selections:
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
    def __init__(self,
                residues: List[str], 
                atoms: Tuple[str],
                segment: Tuple[int] = None,
                *,
                proR: bool = False,
                proS: bool = False) \
                -> None:

        self.plus_minus, split_atoms = self._get_plus_minus_atoms(atoms)
        super().__init__(residues, split_atoms, segment)
        self.proR = proR
        self.proS = proS

    @property
    def selections(self) -> Selections:
        residue_selections = {}
        for residue in self.residues:
            selected_atoms = []
            for atom in self.atoms:
                if atom in constants.SIMPLE_ATOMS:
                    selected_atoms.append([atom])
                elif residue in constants.METHYL_ATOMS.keys() and \
                     (atom == 'HMETHYL' or atom == 'CMETHYL'):
                    if self.proR:
                        methyl_selection = constants.METHYL_ATOMS_PROR
                    elif self.proS:
                        methyl_selection = constants.METHYL_ATOMS_PROS
                    else:
                        methyl_selection = constants.METHYL_ATOMS
                    # first atom is always proton, second is carbon in METHYL_ATOMS dicts
                    index = 0 if atom == 'HMETHYL' else 1
                    selected_atoms.append([methyl_atoms[index] 
                                          for methyl_atoms
                                          in methyl_selection[residue]])
                elif atom == 'H#':
                    # needs a separate string since H is a valid atom (Hn)
                    selected_atoms.append([residue_atom 
                                          for residue_atom
                                          in constants.RESIDUE_ATOMS[residue] 
                                          if residue_atom.startswith('H')])
                elif atom == 'C#':
                    # needs a separate string since C is a valid atom (Co)
                    selected_atoms.append([residue_atom 
                                          for residue_atom
                                          in constants.RESIDUE_ATOMS[residue] 
                                          if residue_atom.startswith('C')])
                else:
                    selected_atoms.append([residue_atom 
                                          for residue_atom
                                          in constants.RESIDUE_ATOMS[residue] 
                                          if residue_atom.startswith(atom)])
            residue_selections[residue] = list(itertools.product(*selected_atoms))
        return residue_selections
    
    def get_correlations(self, entity_shifts: List[NamedTuple], *, label=None) -> Correlations:
        # Intra-residue correlation same as the simple case
        if self.plus_minus == [0]*len(self.atoms):
            return super().get_correlations(entity_shifts)

        # Inter-residue correlation
        else:
            # Group atoms by residue, each item is all assigned shifts for that residue
            residue_shifts = self._get_residue_shifts(entity_shifts)
            shift_table = []
            
            # Filter out residues if sequence numbers for i+- cannot be found
            for ires in residue_shifts:
                window = tuple(ires + index for index in self.plus_minus)
                residues = [residue_shifts[res_num] 
                           for res_num in window
                           if residue_shifts.get(res_num) is not None]
                
                # Get residue number, type and chemical shifts for each atom in group
                group = []
                for i, residue in enumerate(residues):
                    selections = self.selections.get(residue[0].Comp_ID)
                    if selections is not None:
                        # Have to filter it to prevent duplications in branched amino acid residues
                        correlations = [s[i] for s in self.selections[residue[0].Comp_ID]]
                        unique_correlations = sorted(set(correlations), key=correlations.index)

                        selected_atoms = []
                        for correlation in unique_correlations:
                            for atom in residue:
                                if atom.Atom_ID == correlation:
                                    selected_atoms.append((int(atom.Seq_ID), atom.Comp_ID, float(atom.Val)))
                        group.append(selected_atoms)

                # Reformat data and filter by those which have all atoms assigned
                for group_data in itertools.product(*group):
                    try:
                        group_sequence_numbers, group_residue_types, group_shifts = zip(*group_data)
                        if label is not None:
                            sequence_number = group_sequence_numbers[label]
                            residue_type = group_residue_types[label]
                        else:
                            ires_index = self.plus_minus.index(0)
                            sequence_number = group_sequence_numbers[ires_index]
                            residue_type = group_residue_types[ires_index]
                        if len(group_shifts) == len(self.atoms):
                            shift_table.append([sequence_number, residue_type, group_shifts])
                    except (IndexError, ValueError):
                        # not all atoms have assigned shifts or residue(s) filtered out
                        continue
            return shift_table
    
    def _get_plus_minus_atoms(self, atoms: str):
        # Advanced correlation setup
        atoms = [a.upper() for a in atoms]
        plus_minus = [0]*len(atoms)
        split_atoms = []
        for i, spin in enumerate(atoms):
            try:
                # Split at plus or minus sign
                atom, sign, index = re.split(r'(\+|\-)', spin)
                split_atoms.append(atom)
                plus_minus[i] = int(sign+index)
            except ValueError:
                # No plus or minus for atom
                split_atoms.append(spin)
                continue
        # Ensure there is an i residue, not all have +/- indices
        if 0 not in plus_minus:
            abs_min = min([abs(pm) for pm in plus_minus])
            add_min = [pm + abs_min for pm in plus_minus]
            if 0 in add_min:
                plus_minus = add_min
                logging.warn(f"No 'i-residue' found. Adjusted all indicies by +{abs_min}.")
            else:
                plus_minus = [pm - abs_min for pm in plus_minus]
                logging.warn(f"No 'i-residue' found. Adjusted all indicies by -{abs_min}.")
        return tuple(plus_minus), tuple(split_atoms)

class GroupSelector:
    def __init__(self,
                residues: List[str], 
                atoms: Tuple[str],
                segment: Tuple[int] = None,
                *,
                proR: bool = False,
                proS: bool = False) \
                -> None:
        
        atom_combinations = list(itertools.product(*atoms))
        self.selectors = [AdvancedSelector(residues, combination, segment, proR=proR, proS=proS)
                         for combination in atom_combinations]
        self.atoms = self.selectors[0].atoms
    
    def get_correlations(self, entity_shifts: List[NamedTuple], *, label=None) -> Correlations:
        shift_table = []
        for selector in self.selectors:
            shift_table.extend(selector.get_correlations(entity_shifts, label=label))
        return shift_table