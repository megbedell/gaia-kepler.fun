from astropy.utils.data import download_file
from astropy.io import ascii
from astropy.table import QTable, Table
from astropy import units as u
from units import kepler_unit_map

STELLAR_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=q1_q17_dr25_stellar')
ALIAS_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=keplernames')
EXOPLANETS_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=exoplanets')
KOIS_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=q1_q17_dr25_koi')
K2TARGETS_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=k2targets')
K2CAND_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=k2candidates')
                      
def get_table(url=None, cache=True, show_progress=True,
                    table_path=None, select=None, format='csv'):
    """
    Download (and optionally cache) a table from the `NExScI Exoplanet Archive 
                                <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.
                                
    Parameters
    ----------
    url : str (optional)
        Web location of the table to be downloaded.
    cache : bool (optional)
        Cache table to local astropy cache? Default is `True`.
    show_progress : bool (optional)
        Show progress of table download (if no cached copy is
        available). Default is `True`.
    table_path : str (optional)
        Path to a local table file. Default `None` will trigger a
        download of the table from the internet.
    select : str (optional)
        Comma-separated, no spaces string indicating columns to be
        returned. Default `None` will select all default columns
        as set by the Exoplanet Archive API.
    Returns
    -------
    table : `~astropy.table`
        Astropy table of requested data.
    """
    if table_path is None:
        url += '&format='+format
        if not select is None:
            url += '&select='+select
        table_path = download_file(url, cache=cache,
                                   show_progress=show_progress,
                                   timeout=120)
    table = ascii.read(table_path)

    # Assign units to columns where possible
    for col in table.colnames:
        if col in kepler_unit_map:
            if not isinstance(kepler_unit_map[col], u.UnrecognizedUnit): # unit is valid
                table[col].unit = kepler_unit_map[col]

    return table

    
def get_confirmed_planets_table(cache=True, show_progress=True,
                                table_path=None, select=None, format='csv'):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive Confirmed
    Planets table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_exoplanet_columns.html>`_
    """
    exoplanet_table = get_table(url=EXOPLANETS_URL, cache=cache, 
                                    show_progress=show_progress, 
                                    table_path=table_path, select=select, format=format)

    
    # Store column of lowercase names for indexing:
    lowercase_names = [host_name.lower().replace(' ', '') + letter
                       for host_name, letter in
                       zip(exoplanet_table['pl_hostname'].data,
                           exoplanet_table['pl_letter'].data)]
    exoplanet_table['pl_name'] = lowercase_names
    exoplanet_table.add_index('pl_name')
    
    
    return exoplanet_table
    
def get_kois_table(cache=True, show_progress=True,
                            table_path=None, select=None, format='csv'):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive Kepler
    Objects of Interest table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html>`_
    """
    koi_table = get_table(url=KOIS_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select, format=format)
    
    return koi_table
    
def get_keplerstellar_table(cache=True, show_progress=True,
                                table_path=None, select=None, format='csv'):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive Kepler
    Stellar table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html>`_
    """
    keplerstellar_table = get_table(url=STELLAR_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select, format=format)
    
    return keplerstellar_table
    
def get_k2targets_table(cache=True, show_progress=True,
                                table_path=None, select=None, format='csv'):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive K2
    Targets table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_k2_columns.html>`_
    """
    k2targets_table = get_table(url=K2TARGETS_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select, format=format)
    
    return k2targets_table
    
def get_k2candidates_table(cache=True, show_progress=True,
                                table_path=None, select=None, format='csv'):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive K2
    Candidates table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_k2candidates_columns.html>`_
    """
    k2candidates_table = get_table(url=K2CAND_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select, format=format)
    
    return k2candidates_table
    
def get_alias_table(cache=True, show_progress=True,
                                table_path=None, select=None, format='csv'):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive Kepler
    Names table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    A full description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_keplernames_columns.html>`_
    """
    alias_table = get_table(url=ALIAS_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select, format=format)
    
    return alias_table
