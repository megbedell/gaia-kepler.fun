import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
from astroquery.xmatch import XMatch
from astropy.table import QTable, Table, Column, MaskedColumn, join, unique, vstack
import pandas as pd
from tqdm import tqdm
import os
import json
from nasa_tables import *
from units import gaia_unit_map
from plots import skyview_cmd
import pdb

gaia_col_keep = ['angDist', 'kepid', 'epic_number', 'pl_name', 'ra_ep2000', 'dec_ep2000', 'ra', 'dec', 'hip', 
                'tycho2_id', 'solution_id', 'source_id', 'ref_epoch', 'ra_error', 'dec_error', 
                'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 
                'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_primary_flag',
                'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_mag', 
                'phot_variable_flag']

def xmatch_cds_from_csv(file, ra_name='RA', dec_name='Dec', dist=5, cat='src'):
    """
    do a CDS XMatch between targets in file and Gaia; return astropy table
    """
    if cat=='tgas':
        gaia_cat = 'vizier:I/337/tgas' # DR1 TGAS
    elif cat=='src':
        gaia_cat = 'vizier:I/337/gaia' # DR1 src
    elif cat=='dr2':
        gaia_cat = 'vizier:I/345/gaia2' # DR2
    print("starting xmatch between {0} and {1} cat".format(file, cat))
    table = XMatch.query(cat1=open(file), cat2=gaia_cat, max_distance=dist * u.arcsec, 
                         colRA1=ra_name, colDec1=dec_name)
    print("xmatch complete")
    # Assign units to columns where possible
    for col in table.colnames:
        if col not in gaia_col_keep and col not in gaia_unit_map: # don't care about this quantity
            table.remove_column(col)
        if col in gaia_unit_map:
            if not isinstance(gaia_unit_map[col], u.UnrecognizedUnit): # unit is valid
                table[col].unit = gaia_unit_map[col]
    return table
    
def plot_matches(xmatch_table, names_col, all_names, basename, dist, cat):
    matches = np.asarray([np.sum(xmatch_table[names_col] == k) for k in tqdm(all_names)])
    plt.hist(matches)
    plt.ylabel('Sources')
    plt.xlabel('Matches')
    plt.savefig('{b}_{c}_matches_{d}arcsec.png'.format(b=basename, d=dist, c=cat))
    plt.clf()
    print("{0} sources have Gaia {1} matches.".format(np.sum(matches > 0), cat))
    print("{0} sources have >1 matches.".format(np.sum(matches > 1), cat))
    print("{0} sources have no matches.".format(np.sum(matches == 0), cat))    

    plt.hist(xmatch_table['angDist'], bins=np.arange(0.,dist,0.1))
    plt.ylabel('Sources')
    plt.xlabel('Angular Dist (arcsec)')
    plt.savefig('{b}_{c}_angdist_{d}arcsec.png'.format(b=basename, d=dist, c=cat))
    plt.clf()
    
def save_output(table, filename):
    """
    take astropy table and write out
    """
    table.write(filename+'.fits', format='fits', overwrite=True)
    #table_header = t.meta
    
def run_kepler_query(dist, cat, remake_csvs=False, make_plots=True):
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
    lc_table = join(lc_xmatch_table, lc_nasa_table, keys='kepid', table_names=['gaia', 'kic'], 
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
    
    try:
        save_output(lc_table, '../data/kepler_{c}_{d}arcsec'.format(d=dist, c=cat))
    except:
        print("save failed.")
        pdb.set_trace()
    return lc_table 
       
def run_k2_query(dist, cat, remake_csvs=False, make_plots=True):    
    # K2 targets:
    select = 'epic_number,tm_name,k2_campaign_str,k2_type,ra,dec,k2_lcflag,k2_scflag'
    select += ',k2_teff,k2_tefferr1,k2_tefferr2,k2_logg,k2_loggerr1,k2_loggerr2'
    select += ',k2_metfe,k2_metfeerr1,k2_metfeerr2,k2_rad,k2_raderr1,k2_raderr2'
    select += ',k2_mass,k2_masserr1,k2_masserr2'
    #select = None
    k2_nasa_table = get_k2targets_table(select=select)
    epic = k2_nasa_table['epic_number']
     
    k2_files = ['epic_coords0.csv', 'epic_coords1.csv', 'epic_coords2.csv', 'epic_coords3.csv', 
                'epic_coords4.csv', 'epic_coords5.csv', 'epic_coords6.csv', 'epic_coords7.csv']  
    if remake_csvs or not os.path.isfile(k2_files[0]): # make the file
        print('making K2 target files...')
        k2_coords_table = get_k2targets_table(select='epic_number,ra,dec')
        ra, dec = k2_coords_table['ra'], k2_coords_table['dec']
        fs = [open(f, 'w') for f in k2_files]
        for f in fs:
            f.write('epic_number, ra_epic, dec_epic\n')
        for e,(r,d) in tqdm(zip(epic, zip(ra, dec))):
            file_pick = np.random.randint(0,len(k2_files))
            fs[file_pick].write('{0}, {1:.8f}, {2:.8f}\n'.format(e, r, d))
        for f in fs:
            f.close()
        print('files made')
    
    xmatch_tables = []
    for f in k2_files:
        xmatch_tables.append(xmatch_cds_from_csv(f, ra_name='ra_epic', dec_name='dec_epic', dist=dist, cat=cat))
    k2_xmatch_table = vstack(xmatch_tables)
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
    
    try:
        save_output(k2_table, '../data/k2_{c}_{d}arcsec'.format(d=dist, c=cat)) 
    except:
        print("save failed.")
        pdb.set_trace()
    return k2_table   
    
def get_confirmed(dist, cat, remake_csvs=False, make_plots=True):    
    # confirmed planets:
    print('assembling confirmed planets table...')
    #select = 'pl_hostname,pl_letter,pl_discmethod,pl_pnum,pl_orbper,pl_orbsmax,pl_orbeccen'
    #select += ',pl_bmassj,pl_radj,pl_dens,pl_eqt,pl_insol'
    select = None
    conf_nasa_table = get_confirmed_planets_table(select=select)
    confirmed_file = 'confirmed_coords.csv'
    
    if remake_csvs or not os.path.isfile(confirmed_file): # make the file
        print('making confirmed target file...')
        ra, dec = conf_nasa_table['ra'], conf_nasa_table['dec']
        f = open(confirmed_file, 'w')
        f.write('pl_name, ra_nasa, dec_nasa\n')
        for k,(r,d) in tqdm(zip(conf_nasa_table['pl_name'], zip(ra, dec))):
            f.write('{0}, {1:.8f}, {2:.8f}\n'.format(k, r, d))
        f.close()
        print('file made')
    
    conf_xmatch_table = xmatch_cds_from_csv(confirmed_file, ra_name='ra_nasa', dec_name='dec_nasa', dist=dist, cat=cat)
    confirmed_table = join(conf_xmatch_table, conf_nasa_table, keys='pl_name', table_names=['gaia', 'nasa'], 
                    join_type='left')
     
    if make_plots:  # plots 
        print('making confirmed planets plots...')
        plot_matches(confirmed_table, 'pl_name', conf_nasa_table['pl_name'], 
                     'confirmed', dist, cat)
                     
    try:
        save_output(confirmed_table, '../data/confirmed_{c}_{d}arcsec'.format(d=dist, c=cat))
    except:
        print("save failed.")
        pdb.set_trace()
    return confirmed_table
    

    
if __name__ == "__main__":
    dist = 4 # arcsec radius of query
    cat = 'dr2' # 'src' or 'tgas' or 'dr2'
    remake_csvs = False
    make_plots = True
    
    if False:
        print("running queries with dist = {0} arcsec".format(dist))
        kepler_table = run_kepler_query(dist, cat, remake_csvs=remake_csvs, make_plots=make_plots)
        k2_table = run_k2_query(dist, cat, remake_csvs=remake_csvs, make_plots=make_plots)
        confirmed_table = get_confirmed(dist, cat, remake_csvs=remake_csvs, make_plots=make_plots)
    else:
        print("loading pre-queried data with dist = {0} arcsec".format(dist))
        kepler_table = Table.read('../data/kepler_{c}_{d}arcsec.fits'.format(d=dist, c=cat), format='fits')
        k2_table = Table.read('../data/k2_{c}_{d}arcsec.fits'.format(d=dist, c=cat), format='fits')
        confirmed_table = Table.read('../data/confirmed_{c}_{d}arcsec.fits'.format(d=dist, c=cat), format='fits')
            
    dist = 1 # arcsec radius of query
      
    if True: 
        print("running queries with dist = {0} arcsec".format(dist)) 
        kepler_table = run_kepler_query(dist, cat, remake_csvs=remake_csvs, make_plots=make_plots)
        k2_table = run_k2_query(dist, cat, remake_csvs=remake_csvs, make_plots=make_plots)
        confirmed_table = get_confirmed(dist, cat, remake_csvs=remake_csvs, make_plots=make_plots)
    
