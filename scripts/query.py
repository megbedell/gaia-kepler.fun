import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
from astropy.utils.data import download_file
from astroquery.xmatch import XMatch
from astropy.table import Table, Column, MaskedColumn, join, unique
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
K2TARGETS_CSV_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=k2targets')
K2CAND_CSV_URL = ('http://exoplanetarchive.ipac.caltech.edu/cgi-bin/'
                      'nstedAPI/nph-nstedAPI?table=k2candidates')
                      
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
    keplerstellar_table = get_table(url=STELLAR_CSV_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select)
    
    return keplerstellar_table
    
def get_k2targets_table(cache=True, show_progress=True,
                                table_path=None, select=None):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive K2
    Targets table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_k2_columns.html>`_
    """
    k2targets_table = get_table(url=K2TARGETS_CSV_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select)
    
    return k2targets_table
    
def get_k2candidates_table(cache=True, show_progress=True,
                                table_path=None, select=None):
    """
    Download (and optionally cache) the `NExScI Exoplanet Archive K2
    Candidates table <http://exoplanetarchive.ipac.caltech.edu/index.html>`_.

    The Exoplanet Archive table returns lots of columns of data. A full
    description of the columns can be found `here
    <https://exoplanetarchive.ipac.caltech.edu/docs/API_k2candidates_columns.html>`_
    """
    k2candidates_table = get_table(url=K2CAND_CSV_URL, cache=cache, 
                                        show_progress=show_progress, 
                                        table_path=table_path, select=select)
    
    return k2candidates_table
    
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
    
def save_plot_data(table, filename, subsample_size=4999):
    df = table.to_pandas()
    N_tot = len(df)
    subsample_ind = np.random.permutation(np.arange(N_tot))[:subsample_size] # random values, no duplication
    subsample = pd.DataFrame(df.iloc[list(subsample_ind)])
    try:
        gmk = subsample['phot_g_mean_mag'] - subsample['kepmag']
    except:
        gmk = subsample['phot_g_mean_mag'] - subsample['k2_kepmag']
    subsample['gaiamag_minus_kepmag'] = gmk
    dist = 1.e3/subsample['parallax'] # distance in pc
    subsample['dist'] = dist
    abs_gmag = subsample['phot_g_mean_mag'] - 5.*(np.log10(dist) - 1.)
    subsample['abs_gmag'] = abs_gmag
    subsample.to_csv(filename+'_subsample.csv', mode='w')
    
if __name__ == "__main__":
    dist = 1 # arcsec radius of query
    cat = 'tgas' # 'src' or 'tgas'
    remake_csvs = False
    make_plots = False
    
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
        
    # flag KOIs and confirmed hosts:
    full_alias_table = get_alias_table()
    alias_table = unique(full_alias_table, keys='kepid')
    alias_table_to_join = alias_table['kepid', 'kepoi_name']
    lc_table = join(lc_table, alias_table_to_join, keys='kepid', join_type='left')
    lc_table['planet?'] = 'none'
    lc_table['planet?'][~lc_table['kepoi_name'].mask] = 'candidate'
    # TODO: flag confirmed
    
    save_output(lc_table, '../data/lc_{c}_{d}arcsec'.format(d=dist, c=cat)) 
    save_plot_data(lc_table, '../data/plot_lc_{c}_{d}arcsec'.format(d=dist, c=cat))    
       
    
    # K2 targets:
    k2_file = 'epic_coords.csv'  
    #select = 'epic_number,tm_name,k2_campaign_str,k2_type,ra,dec,k2_lcflag,k2_scflag'
    #select += ',k2_teff,k2_tefferr1,k2_tefferr2,k2_logg,k2_loggerr1,k2_loggerr2'
    #select += ',k2_metfe,k2_metfeerr1,k2_metfeerr2,k2_rad,k2_raderr1,k2_raderr2'
    #select += ',k2_mass,k2_masserr1,k2_masserr2'
    select = None
    k2_nasa_table = get_k2targets_table(select=select)
    epic = k2_nasa_table['epic_number']
     
    if remake_csvs: # make the file
        print('making K2 target file...')
        k2_coords_table = get_k2targets_table(select='epic_number,ra,dec')
        ra, dec = k2_coords_table['ra'], k2_coords_table['dec']
        f = open(k2_file, 'w')
        f.write('epic_number, RA, Dec\n')
        for e,(r,d) in tqdm(zip(epic, zip(ra, dec))):
            f.write('{0}, {1:.8f}, {2:.8f}\n'.format(e, r, d))
        f.close()
        print('file made')
    

    k2_xmatch_table = xmatch_cds_from_csv(k2_file, ra_name='RA', dec_name='Dec', dist=dist, cat=cat)
    k2_table = join(k2_xmatch_table, k2_nasa_table, keys='epic_number', table_names=['gaia', 'nasa'], 
                    join_type='left')
    
    if make_plots:  # plots 
        print('making K2 plots...')
        k2_matches = [np.sum(k2_table['epic_number'] == e) for e in tqdm(epic)]
        plt.hist(k2_matches)
        plt.ylabel('Sources')
        plt.xlabel('Matches')
        plt.savefig('k2_{c}_matches_{d}arcsec.png'.format(d=dist, c=cat))
        plt.clf()
    
        plt.hist(k2_table['angDist'], bins=np.arange(0.,dist,0.1))
        plt.ylabel('Sources')
        plt.xlabel('Angular Dist (arcsec)')
        plt.savefig('k2_{c}_angdist_{d}arcsec.png'.format(d=dist, c=cat))
        plt.clf()
        
    # flag KOIs and confirmed hosts:
    full_k2cand_table = get_k2candidates_table()
    k2cand_table = unique(full_k2cand_table, keys='epic_name')
    k2cand_table_to_join = k2cand_table['k2c_disp','k2c_note']
    epic_numbers = [int(n.split(" ")[1]) for n in k2cand_table['epic_name']]
    k2cand_table_to_join['epic_number'] = epic_numbers
    k2_table = join(k2_table, k2cand_table_to_join, keys='epic_number', join_type='left')
    
    save_output(k2_table, '../data/k2_{c}_{d}arcsec'.format(d=dist, c=cat))    
    save_plot_data(k2_table, '../data/plot_k2_{c}_{d}arcsec'.format(d=dist, c=cat))    
    
    
    # confirmed planets:
    print('assembling confirmed planets table...')
    select = 'pl_hostname,pl_letter,pl_discmethod,pl_pnum,pl_orbper,pl_orbsmax,pl_orbeccen'
    select += ',pl_bmassj,pl_radj,pl_dens,pl_eqt,pl_insol'
    confirmed_nasa_table = get_confirmed_planets_table(select=select)
    
    star_names = [n.split(" ")[0] for n in alias_table['alt_name']] # HACK - might be losing some?
    alias_table['pl_hostname'] = star_names
    alias_table.remove_column('alt_name')
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
    
    
    
    
        
