from collections import namedtuple
import logging
import shutil
from typing import List, NamedTuple
import webbrowser

import requests
from requests.api import head
from requests.exceptions import HTTPError

from nightshift import constants

class BMRBEntity:
    def __init__(self, entry_id: str, name: str, sequence: str, shifts: List[NamedTuple]) -> None:
        self.entry_id = entry_id
        self.name = name
        self.sequence = sequence.replace('\n', '')
        self.shifts = shifts
    
    @property
    def molecular_weight(self) -> float:
        '''Used to estimate the linewidth for contour plots. Not all entries have an
        Entity.Formula_weight tag, but all have a sequence to calculate it from.
        Weights from: https://github.com/biopython/biopython/blob/master/Bio/Data/IUPACData.py'''

        # account for peptide bond formation
        water_weight = 18.0153 * (len(self.sequence) - 1)
        return sum(constants.AMINO_ACID_WEIGHTS[aa] for aa in self.sequence) - water_weight

def get_bmrb_shifts(entry_id: str) -> BMRBEntity:
    '''Use BMRB API to get entry number's entity names and assigned chemical shifts
    Returns None if entry cannot be found or rate limited.'''
    
    bmrb_entry_url = get_bmrb_url('current', 'entry', entry_id)
    
    # Cannot requests tags and loops at the same time, so two requests
    # We could get multiple loops and tags per request. Just using one currently.
    tags_payload = {'tag': ['Entity.name', 'Entity.Polymer_seq_one_letter_code']}
    loops_payload = {'loop': 'Atom_chem_shift',}
    headers = {'Application': 'nightshift'}

    # Parses JSON returned by BMRB into a dictionary
    with requests.get(bmrb_entry_url, params=tags_payload, headers=headers) as tag_request:
        if _got_bad_code(tag_request):
            return
        names = tag_request.json()[entry_id]['Entity.name']
        sequences = tag_request.json()[entry_id]['Entity.Polymer_seq_one_letter_code']

    # Select which entity to plot if more than one
    entity_selection = 0
    if len(names) > 1:
        entity_selection = choose_entity(names)

    # Use BMRB API to get entry number's assigned chemical shifts
    with requests.get(bmrb_entry_url, loops_payload) as loop_request:
        if _got_bad_code(loop_request):
            return
        loop_json = loop_request.json()
        nmrstar = loop_json[entry_id]['Atom_chem_shift'][entity_selection]['data']
        tags = loop_json[entry_id]['Atom_chem_shift'][entity_selection]['tags']

        AssignedChemicalShift = namedtuple('AssignedChemicalShift', tags)
        # List of named tuples indexed by NMRStar field names
        shifts = [AssignedChemicalShift(*row) for row in nmrstar]
    entity = BMRBEntity(entry_id, names[entity_selection], sequences[entity_selection], shifts)
    return entity

def get_bmrb_url(version: str, *args: str) -> str:
    '''Convenience function for making different API calls more easily.'''
    return f'http://api.bmrb.io/{version}/' + '/'.join(args)

def _got_bad_code(response: requests.Response) -> bool:
    try:
        response.raise_for_status()
    except HTTPError as err:
        if response.status_code == 403:
            print('Rate limited by server. Wait 10 seconds to try again.')
            logging.debug('Rate limited ' + err)
        else:
            message = response.json()['error']
            print(message)
            logging.debug(err)
        return True
    return False

def choose_entity(names: list) -> int:
    while True:
        try:
            print('Found more than one assigned chemical shift list, please select:')
            # Print list of names received from names_url request number matches Entity_ID tag
            for i, name in enumerate(names, start=1):
                print(f'({i}) {name}')
            entity_selection = int(input('> ')) - 1 # Entity_ID is 1-indexed, python is 0-indexed
            if 0 <= entity_selection < len(names):
                return entity_selection
            else:
                print(f'Invalid selection please choose a number <= {len(names)}.')
        except ValueError:
            print(f'Invalid selection, please input a number <= {len(names)}.')

def bmrb_instant_search(terms: str) -> None:
    bmrb_search_url = get_bmrb_url('current', 'instant')
    search_payload = {'term': terms, 'database': 'macromolecules'}
    headers = {'Application': 'nightshift'}

    terminal_width = shutil.get_terminal_size()[0]
    padding = 8
    with requests.get(bmrb_search_url, params=search_payload, headers=headers) as search_request:
        if _got_bad_code(search_request):
            return
        all_results = search_request.json()
        if all_results:
            for result in all_results:
                # Only print results with assigned chemical shifts
                data = set(r.get('type') for r in result['data_types'])
                if 'assigned_chemical_shifts' in data:
                    line = f"{result['value']}\t{result['label']}"
                    if len(line) + padding > terminal_width:
                        line = line[:terminal_width - padding] + '...'
                    print(line)
        else:
            print(f"No entries with assigned chemical shifts for '{terms}'")

def open_bmrb_page(entry_id) -> None:
    # First check to see if the entry number is valid, then open in webbrowser
    bmrb_entry_url = get_bmrb_url('current', 'entry', entry_id)
    headers = {'Application': 'nightshift'}
    with requests.get(bmrb_entry_url, headers=headers) as webpage_request:
        if not _got_bad_code(webpage_request):
            webbrowser.open(f'https://bmrb.io/data_library/summary/index.php?bmrbId={entry_id}')