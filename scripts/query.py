import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
from astropy.utils.data import download_file
from astroquery.xmatch import XMatch
from astropy.table import Table, Column, MaskedColumn, join
import pandas as pd
from tqdm import tqdm

STELLAR_CSV_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=q1_q17_dr25_stellar')
ALIAS_CSV_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=keplernames')
EXOPLANETS_CSV_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=exoplanets')
KOIS_CSV_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=q1_q17_dr25_koi')

                      
def get_table(url=None, cache=True, show_progress=True,
                    table_path=None, select=None):
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
        if not select is None:
            url += '&select='+select
        table_path = download_file(url, cache=cache,
                                   show_progress=show_progress,
                                   timeout=120)
    table = ascii.read(table_path)

    return table

    
def get_confirmed_planets_table(cache=True, show_progress=True,
                                table_path=None, select=None):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive Confirmed
    Planets table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_exoplanet_columns.html>`_
    """
    exoplanet_table = get_table(url=EXOPLANETS_CSV_URL, cache=cache, 
                                    show_progress=show_progress, 
                                    table_path=table_path, select=select)

    
    # Store column of lowercase names for indexing:
    lowercase_names = [host_name.lower().replace(' ', '') + letter
                       for host_name, letter in
                       zip(exoplanet_table['pl_hostname'].data,
                           exoplanet_table['pl_letter'].data)]
    exoplanet_table['pl_name'] = lowercase_names
    exoplanet_table.add_index('pl_name')
    
    
    return exoplanet_table
    
def get_kois_table(cache=True, show_progress=True,
                            table_path=None, select=None):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive Kepler
    Objects of Interest table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html>`_
    """
    select = None
    koi_table = get_table(url=KOIS_CSV_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select)
    
    return koi_table
    
def get_keplerstellar_table(cache=True, show_progress=True,
                                table_path=None, select=None):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive Kepler
    Stellar table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html>`_
    """
    select = None
    keplerstellar_table = get_table(url=STELLAR_CSV_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select)
    
    return keplerstellar_table
    
def get_alias_table(cache=True, show_progress=True,
                                table_path=None, select=None):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive Kepler
    Names table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    A full description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_keplernames_columns.html>`_
    """
    alias_table = get_table(url=ALIAS_CSV_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select)
    
    return alias_table

def xmatch_cds_from_csv(file, ra_name='RA', dec_name='Dec', dist=5, cat='src'):
    """
    do a CDS XMatch between targets in file and Gaia; return astropy table
    """
    if cat=='tgas':
        gaia_cat = 'vizier:I/337/tgas' # DR1 TGAS
    elif cat=='src':
        gaia_cat = 'vizier:I/337/gaia'
    table = XMatch.query(cat1=open(file), cat2=gaia_cat, max_distance=dist * u.arcsec, colRA1=ra_name, colDec1=dec_name)
    return table
    
def save_output(table, filename):
    # take astropy table and write out
    table.write(filename+'.fits', format='fits', overwrite=True)
    df = table.to_pandas()
    df.to_hdf(filename+'.h5', 'data', mode='w')
    
if __name__ == "__main__":
    dist = 1 # arcsec radius of query
    cat = 'src' # 'src' or 'tgas'
    remake_csvs = False
    make_plots = True
    
    # OG Kepler long-cadence targets:
    lc_file = 'kic_lc_coords.csv'  
    lc_nasa_table = get_keplerstellar_table()
    kic = lc_nasa_table['kepid']
     
    if remake_csvs: # make the file
        print('making LC target file...')
        ra, dec = lc_nasa_table['ra'], lc_nasa_table['dec']
        f = open(lc_file, 'w')
        f.write('kepid, RA, Dec\n')
        for k,(r,d) in tqdm(zip(kic, zip(ra, dec))):
            f.write('{0}, {1:.8f}, {2:.8f}\n'.format(k, r, d))
        f.close()
        print('file made')
    

    lc_xmatch_table = xmatch_cds_from_csv(lc_file, ra_name='RA', dec_name='Dec', dist=dist, cat=cat)
    lc_table = join(lc_xmatch_table, lc_nasa_table, keys='kepid', table_names=['gaia', 'nasa'], 
                    join_type='left')
    save_output(lc_table, '../data/lc_{c}_{d}arcsec'.format(d=dist, c=cat))    
    
    if make_plots:  # plots 
        print('making LC plots...')
        lc_matches = [np.sum(lc_table['kepid'] == k) for k in tqdm(kic)]
        plt.hist(lc_matches)
        plt.ylabel('Sources')
        plt.xlabel('Matches')
        plt.savefig('lc_{c}_matches_{d}arcsec.png'.format(d=dist, c=cat))
        plt.clf()
    
        plt.hist(lc_table['angDist'], bins=np.arange(0.,dist,0.1))
        plt.ylabel('Sources')
        plt.xlabel('Angular Dist (arcsec)')
        plt.savefig('lc_{c}_angdist_{d}arcsec.png'.format(d=dist, c=cat))
        plt.clf()
    
    
    # OG Kepler short-cadence targets:
    print('assembling SC table...')
    sc_file = 'kepler_sc_kicnames.txt'
    sc_names = np.genfromtxt(sc_file, delimiter=',', unpack=True, dtype='int64')
    sc_table = lc_table.copy()
    delete = []
    for i in tqdm(range(len(sc_table))):
        if not sc_table['kepid'][i] in sc_names:
            delete.append(i)
    sc_table.remove_rows(delete)
    save_output(sc_table, '../data/sc_{c}_{d}arcsec'.format(d=dist, c=cat))
    
    
    if make_plots:
        print('making SC plots...')
        sc_matches = [np.sum(sc_table['kepid'] == k) for k in tqdm(sc_names)]
        plt.hist(sc_matches)
        plt.ylabel('Sources')
        plt.xlabel('Matches')
        plt.savefig('sc_{c}_matches_{d}arcsec.png'.format(d=dist, c=cat))
        plt.clf()
    
    
        plt.hist(sc_table['angDist'], bins=np.arange(0.,dist,0.1))
        plt.ylabel('Sources')
        plt.xlabel('Angular Dist (arcsec)')
        plt.savefig('sc_{c}_angdist_{d}arcsec.png'.format(d=dist, c=cat))
        plt.clf()
        
    # OG Kepler KOIs:
    #kois_table = NasaExoplanetArchive.get_kois_table()
    #alias_table = get_alias_table(select='kepid,kepoi_name')
    
    
    # OG Kepler confirmed planets:
    print('assembling Kepler confirmed planets table...')
    select = 'pl_hostname,pl_letter,pl_discmethod,pl_pnum,pl_orbper,pl_orbsmax,pl_orbeccen'
    select += ',pl_bmassj,pl_radj,pl_dens,pl_eqt,pl_insol'
    confirmed_nasa_table = get_confirmed_planets_table(select=select)
    
    alias_table = get_alias_table(select='kepid,alt_name')
    star_names = [n.split(" ")[0] for n in alias_table['alt_name']] # HACK - might be losing some?
    # TODO: test/debug on KOI-351
    alias_table['pl_hostname'] = star_names
    confirmed_nasa_table_with_kepid = join(confirmed_nasa_table, alias_table, keys='pl_hostname',
                                join_type='inner') # add KIC numbers, remove non-Kepler hosts
    confirmed_table = join(confirmed_nasa_table_with_kepid, lc_table, keys='kepid', 
                            join_type='left') # add Gaia results
    save_output(confirmed_table, '../data/confirmed_{c}_{d}arcsec'.format(d=dist, c=cat))
    
    if make_plots:  # plots 
        print('making confirmed planets plots...')
        name = confirmed_nasa_table['pl_name']
        confirmed_matches = [np.sum(confirmed_table['pl_name'] == n) for n in tqdm(name)]
        plt.hist(confirmed_matches)
        plt.ylabel('Sources')
        plt.xlabel('Matches')
        plt.savefig('confirmed_{c}_matches_{d}arcsec.png'.format(d=dist, c=cat))
        plt.clf()
    
        plt.hist(confirmed_table['angDist'], bins=np.arange(0.,dist,0.1))
        plt.ylabel('Sources')
        plt.xlabel('Angular Dist (arcsec)')
        plt.savefig('confirmed_{c}_angdist_{d}arcsec.png'.format(d=dist, c=cat))
        plt.clf()
    

    
    # Kepler FOV:
    
    
    
        
