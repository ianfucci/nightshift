from collections import namedtuple
import logging
from typing import List, NamedTuple, Tuple
import webbrowser

import requests
from requests.exceptions import HTTPError

def get_shifts(entry_id: str) -> Tuple[List, List[List[NamedTuple]]]:
    '''Use BMRB API to get entry number's entity names and assigned chemical shifts.
    Returns a list of entity names and a list of their shifts as a named tuple.
    Returns None if entry cannot be found or rate limited.
    '''
    
    bmrb_entry_url = get_bmrb_url('v2', 'entry', entry_id)
    
    # Cannot request tags and loops at the same time, so two requests
    # We could get multiple loops and tags per request. Just using one currently.
    tags_payload = {'tag': 'Entity.name'}
    loops_payload = {'loop': 'Atom_chem_shift'}
    headers = {'Application': 'nightshift'}

    with requests.get(bmrb_entry_url, params=tags_payload, headers=headers) as tag_request:
        if _got_bad_code(tag_request):
            return
        names = tag_request.json()[entry_id]['Entity.name']

    # Use BMRB API to get entry number's assigned chemical shifts
    with requests.get(bmrb_entry_url, loops_payload) as loop_request:
        if _got_bad_code(loop_request):
            return
        loop_json = loop_request.json()
        # import pdb; pdb.set_trace()
        shift_table = loop_json[entry_id]['Atom_chem_shift']
        data = [entity['data'] for entity in shift_table]
        tags = [entity['tags'] for entity in shift_table]

        # Named tuples indexed by NMRStar field names
        AssignedChemicalShift = namedtuple('AssignedChemicalShift', tags[0])
        shifts = []
        for entity in data:
            entity_shifts = []
            for row in entity:
                entity_shifts.append(AssignedChemicalShift(*row))
            shifts.append(entity_shifts)

    BMRBEntryData = namedtuple('BMRBEntryData', ['names', 'shifts'])
    return BMRBEntryData(names, shifts)

def get_bmrb_url(version: str, *args: str) -> str:
    '''Convenience function for reaching different BMRB API endpoints.
    Additional parameters should be passed as a dictionary to requests
    '''
    return f'http://api.bmrb.io/{version}/' + '/'.join(args)

def _got_bad_code(response: requests.Response) -> bool:
    '''Checks for a 403 response which means >50 requests per second or a 
    404 response which means the entry could not be found.
    '''
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

def instant_search(terms: str) -> List[Tuple[str, str]]:
    '''Uses the BMRB instant search to find the search terms. Can be keyword
    author, PDB code, etc. Returns a list of tuples where the entry ID is the
    first item and the title of the entry is the second.
    '''
    bmrb_search_url = get_bmrb_url('v2', 'instant')
    search_payload = {'term': terms, 'database': 'macromolecules'}
    headers = {'Application': 'nightshift'}

    with requests.get(bmrb_search_url, params=search_payload, headers=headers) as search_request:
        if _got_bad_code(search_request):
            return
        all_results = search_request.json()
        results_with_shifts = []
        if all_results:
            for result in all_results:
                # Only save results with assigned chemical shifts
                data = set(field.get('type') for field in result['data_types'])
                if 'assigned_chemical_shifts' in data:
                    results_with_shifts.append((result['value'], result['label']))
        return results_with_shifts


def open_website(entry_id) -> None:
    '''Opens the BMRB webpage for an entry in the user's default browser.
    '''
    # First check to see if the entry number is valid with a request
    bmrb_entry_url = get_bmrb_url('v2', 'entry', entry_id)
    headers = {'Application': 'nightshift'}
    with requests.get(bmrb_entry_url, headers=headers) as webpage_request:
        if not _got_bad_code(webpage_request):
            webbrowser.open(f'https://bmrb.io/data_library/summary/index.php?bmrbId={entry_id}')