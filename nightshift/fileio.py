from collections import namedtuple
from typing import List, NamedTuple, TextIO

def read_nef_converted_nmrstar(nmrstar: TextIO) -> NamedTuple:
    line = nmrstar.readline()
    while not line.strip().startswith('save_nef_chemical_shift_list'):
        line = nmrstar.readline()
        if not line:
            raise EOFError('NEF format not found')
    
    # enter nef_chemical_shift_list loop
    while not line.strip().startswith('loop_'):
        line = nmrstar.readline()
    line = nmrstar.readline()

    # get fieldnames
    fields = []
    while line.strip().startswith('_nef_chemical_shift.'):
        fields.append(line.strip().split('.')[-1])
        line = nmrstar.readline()

    NEFChemicalShift = namedtuple('NEFChemicalShift', fields)
    shifts = []
    for shift_line in nmrstar:
        if shift_line.strip() == 'stop_':
            break
        if shift_line.strip():
            shifts.append(NEFChemicalShift(*shift_line.strip().split()))
    return make_nmrstar_like(shifts)

def read_shiftx2_output_nmrstar(nmrstar: TextIO) -> NamedTuple:
    # starts at data, no header
    shifts =[]
    ShiftX2ChemicalShift = namedtuple('ShiftX2ChemicalShift', ['ID','Seq_ID','Comp_ID','Atom_ID','Val','Val_err','Entity_ID'])
    for shift_line in nmrstar:
        if shift_line.strip():
            shifts.append(ShiftX2ChemicalShift(*shift_line.strip().split()))
    return make_nmrstar_like(shifts)

def make_nmrstar_like(shifts: List[NamedTuple]) -> List[NamedTuple]:
    # only have oen type so far don't check
    TranslatedChemicalShifts = namedtuple('TranslatedChemicalShifts', ['Seq_ID', 'Comp_ID', 'Atom_ID', 'Val'])
    translated_shifts = []
    for shift in shifts:
        if shift.__class__.__name__ == 'NEFChemicalShift':
            translated = TranslatedChemicalShifts(shift.sequence_code, shift.residue_name, shift.atom_name, shift.value)
            translated_shifts.append(translated)
        if shift.__class__.__name__ == 'ShiftX2ChemicalShift':
            translated = TranslatedChemicalShifts(shift.Seq_ID, shift.Comp_ID, shift.Atom_ID, shift.Val)
            translated_shifts.append(translated)
    return translated_shifts