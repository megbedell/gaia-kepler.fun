import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
from astroquery.xmatch import XMatch
from astropy.table import Table, Column, MaskedColumn, join, unique
import pandas as pd
from tqdm import tqdm
import os
from nasa_tables import *
from plots import skyview_cmd

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
    
def plot_matches(xmatch_table, names_col, all_names, basename, dist, cat):
    matches = [np.sum(xmatch_table[names_col] == k) for k in tqdm(all_names)]
    plt.hist(matches)
    plt.ylabel('Sources')
    plt.xlabel('Matches')
    plt.savefig('{b}_{c}_matches_{d}arcsec.png'.format(b=basename, d=dist, c=cat))
    plt.clf()

    plt.hist(xmatch_table['angDist'], bins=np.arange(0.,dist,0.1))
    plt.ylabel('Sources')
    plt.xlabel('Angular Dist (arcsec)')
    plt.savefig('{b}_{c}_angdist_{d}arcsec.png'.format(b=basename, d=dist, c=cat))
    plt.clf()
    
def save_output(table, filename):
    # take astropy table and write out
    table.write(filename+'.fits', format='fits', overwrite=True)

    
if __name__ == "__main__":
    dist = 1 # arcsec radius of query
    cat = 'tgas' # 'src' or 'tgas'
    remake_csvs = False
    make_plots = False
    
    # OG Kepler long-cadence targets:
    lc_file = 'kic_lc_coords.csv'  
    select = 'kepid,tm_designation,ra,dec,kepmag'
    select += ',teff,teff_err1,teff_err2,teff_prov,logg,logg_err1,logg_err2,logg_prov'
    select += ',feh,feh_err1,feh_err2,feh_prov,radius,radius_err1,radius_err2'
    select += ',mass,mass_err1,mass_err2,prov_sec,nconfp,nkoi,ntce,jmag,hmag,kmag'
    #select=None
    lc_nasa_table = get_keplerstellar_table(select=select)
    kic = lc_nasa_table['kepid']
     
    if remake_csvs or not os.path.isfile(lc_file): # make the file
        print('making LC target file...')
        ra, dec = lc_nasa_table['ra'], lc_nasa_table['dec']
        f = open(lc_file, 'w')
        f.write('kepid, ra_kic, dec_kic\n')
        for k,(r,d) in tqdm(zip(kic, zip(ra, dec))):
            f.write('{0}, {1:.8f}, {2:.8f}\n'.format(k, r, d))
        f.close()
        print('file made')
    
    
    lc_xmatch_table = xmatch_cds_from_csv(lc_file, ra_name='ra_kic', dec_name='dec_kic', dist=dist, cat=cat)
    lc_nasa_table.remove_columns(['ra', 'dec'])
    lc_table = join(lc_xmatch_table, lc_nasa_table, keys='kepid', table_names=['gaia', 'nasa'], 
                    join_type='left')
    
    if make_plots:  # plots 
        print('making LC plots...')
        plot_matches(lc_table, 'kepid', kic, 'lc', dist, cat)
        
    # flag KOIs and confirmed hosts:
    full_alias_table = get_alias_table()
    alias_table = unique(full_alias_table, keys='kepid')
    alias_table_to_join = alias_table['kepid', 'kepoi_name']
    lc_table = join(lc_table, alias_table_to_join, keys='kepid', join_type='left')
    lc_table['planet?'] = 'none'
    lc_table['planet?'][lc_table['nkoi'] > 0] = 'cand'
    lc_table['planet?'][lc_table['nconfp'] > 0] = 'conf'
    
    save_output(lc_table, '../data/lc_{c}_{d}arcsec'.format(d=dist, c=cat)) 
       
    
    # K2 targets:
    k2_file = 'epic_coords.csv'  
    select = 'epic_number,tm_name,k2_campaign_str,k2_type,ra,dec,k2_lcflag,k2_scflag'
    select += ',k2_teff,k2_tefferr1,k2_tefferr2,k2_logg,k2_loggerr1,k2_loggerr2'
    select += ',k2_metfe,k2_metfeerr1,k2_metfeerr2,k2_rad,k2_raderr1,k2_raderr2'
    select += ',k2_mass,k2_masserr1,k2_masserr2'
    #select = None
    k2_nasa_table = get_k2targets_table(select=select)
    epic = k2_nasa_table['epic_number']
     
    if remake_csvs or not os.path.isfile(k2_file): # make the file
        print('making K2 target file...')
        k2_coords_table = get_k2targets_table(select='epic_number,ra,dec')
        ra, dec = k2_coords_table['ra'], k2_coords_table['dec']
        f = open(k2_file, 'w')
        f.write('epic_number, ra_epic, dec_epic\n')
        for e,(r,d) in tqdm(zip(epic, zip(ra, dec))):
            f.write('{0}, {1:.8f}, {2:.8f}\n'.format(e, r, d))
        f.close()
        print('file made')
    

    k2_xmatch_table = xmatch_cds_from_csv(k2_file, ra_name='ra_epic', dec_name='dec_epic', dist=dist, cat=cat)
    k2_table = join(k2_xmatch_table, k2_nasa_table, keys='epic_number', table_names=['gaia', 'nasa'], 
                    join_type='left')
    
    if make_plots:  # plots 
        print('making K2 plots...')
        plot_matches(k2_table, 'epic_number', epic, 'k2', dist, cat)
        
    # flag KOIs and confirmed hosts:
    full_k2cand_table = get_k2candidates_table()
    k2cand_table = unique(full_k2cand_table, keys='epic_name')
    k2cand_table_to_join = k2cand_table['k2c_disp','k2c_note']
    epic_numbers = [int(n.split(" ")[1]) for n in k2cand_table['epic_name']]
    k2cand_table_to_join['epic_number'] = epic_numbers
    k2_table = join(k2_table, k2cand_table_to_join, keys='epic_number', join_type='left')
    
    save_output(k2_table, '../data/k2_{c}_{d}arcsec'.format(d=dist, c=cat))    
    
    
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
        plot_matches(confirmed_nasa_table, 'pl_name', confirmed_nasa_table['pl_name'], 
                     'confirmed', dist, cat)   
                         
