from astropy.table import QTable, Table, Column, MaskedColumn, join, unique, vstack
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
import astropy.coordinates as coord
from astropy.time import Time
from units import gaia_unit_map, kepler_unit_map
from nasa_tables import *

def clean_gaia_table(tbl, kepler=False, k2=False,exoplanets=False):
    """
    Add units, delete some columns
    """
    cols_to_delete = ['solution_id', 'random_index', 'astrometric_n_obs_al', 'astrometric_n_obs_ac', 
                      'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al',
                      'astrometric_params_solved', 'astrometric_weight_al', 'astrometric_pseudo_colour',
                      'astrometric_pseudo_colour_error', 'mean_varpi_factor_al', 'astrometric_matched_observations',
                      'visibility_periods_used', 'astrometric_sigma5d_max',  
                      'frame_rotator_object_type', 'matched_observations', 'phot_g_n_obs', 
                      'phot_g_mean_flux_over_error', 'phot_bp_n_obs', 'phot_bp_mean_flux_over_error',
                      'phot_rp_n_obs', 'phot_rp_mean_flux_over_error', 'phot_bp_rp_excess_factor', 
                      'phot_proc_mode', 'rv_nb_transits', 'rv_template_teff', 'rv_template_logg', 
                      'rv_template_fe_h', 'priam_flags', 'flame_flags', 'datalink_url', 'epoch_photometry_url'
                    ]
    tbl.remove_columns(cols_to_delete)
    if kepler:
        cols_to_delete = ['kepler_oid', 'angdist']
        tbl.remove_columns(cols_to_delete)
    if k2:
        cols_to_delete = ['k2_oid', 'angdist']
        tbl.remove_columns(cols_to_delete)
    if exoplanets:
        cols_to_delete = ['exoplanets_oid', 'angdist']
        tbl.remove_columns(cols_to_delete)
        for i,p in enumerate(tbl['pl_name']):
            tbl['pl_name'][i] = p.strip() # remove whitespace
    for col in tbl.colnames:
        if col in gaia_unit_map:
            if not isinstance(gaia_unit_map[col], u.UnrecognizedUnit): # unit is valid
                tbl[col].unit = gaia_unit_map[col]    
    tbl.rename_column('ref_epoch', 'gaia_ref_epoch')     
    return tbl
    
def clean_dist_table(tbl):
    """
    Add units, delete some columns
    """
    tbl = unique(tbl)
    tbl.rename_column('r_len', 'r_length_prior')
    tbl.rename_column('result_flag', 'r_result_flag')
    tbl.rename_column('modality_flag', 'r_modality_flag')
    for col in ['r_est', 'r_lo', 'r_hi', 'r_length_prior']:
        tbl[col].unit = u.parsec
    return tbl
    
def clean_kepler_table(tbl):
    """
    Add units
    """
    for col in tbl.colnames:
        if col in kepler_unit_map:
            if not isinstance(kepler_unit_map[col], u.UnrecognizedUnit): # unit is valid
                tbl[col].unit = kepler_unit_map[col]  
    return tbl
    
def make_full_tables(data_dir='../data/',kepler=False,k2=False,exoplanets=False):
    """
    Load up NASA table, gaia x-match, and distance x-match
    Combine and write out final tables
    """
    if kepler:
        gaia_matches_file = data_dir+'kepler_5arcsec_gaia.fits'
        dist_table_file = data_dir+'kepler_5arcsec_dist.fits'
        outfile_4arcsec = 'kepler_dr2_4arcsec.fits'
        outfile_1arcsec = 'kepler_dr2_1arcsec.fits'
        
        select = 'kepid,tm_designation,kepmag'
        select += ',teff,teff_err1,teff_err2,teff_prov,logg,logg_err1,logg_err2,logg_prov'
        select += ',feh,feh_err1,feh_err2,feh_prov,radius,radius_err1,radius_err2'
        select += ',mass,mass_err1,mass_err2,prov_sec,nconfp,nkoi,ntce,jmag,hmag,kmag'
        nasa_table = get_keplerstellar_table(select=select, cache=False)
        nasa_table = clean_kepler_table(nasa_table)
        
        # make custom planet host status column
        nasa_table['planet?'] = 'none'
        nasa_table['planet?'][nasa_table['nkoi'] > 0] = 'cand'
        nasa_table['planet?'][nasa_table['nconfp'] > 0] = 'conf'
        
        nasa_table_key = 'kepid' # for join to nasa_table
        ra_key, dec_key = 'ra_kic', 'dec_kic' # for coordinate transforms
        ang_dist_key = 'kepler_gaia_ang_dist' # name for gaia - nasa angular distance
              
    elif k2:
        gaia_matches_file = data_dir+'k2_20arcsec_gaia.fits'
        dist_table_file = data_dir+'k2_20arcsec_dist.fits'
        outfile_20arcsec = 'k2_dr2_20arcsec.fits'
        outfile_4arcsec = 'k2_dr2_4arcsec.fits'
        outfile_1arcsec = 'k2_dr2_1arcsec.fits'
        
        select = 'epic_number,tm_name,k2_campaign_str,k2_type,k2_lcflag,k2_scflag'
        select += ',k2_teff,k2_tefferr1,k2_tefferr2,k2_logg,k2_loggerr1,k2_loggerr2'
        select += ',k2_metfe,k2_metfeerr1,k2_metfeerr2,k2_rad,k2_raderr1,k2_raderr2'
        select += ',k2_mass,k2_masserr1,k2_masserr2,k2_kepmag,k2_kepmagerr,k2_kepmagflag'
        nasa_table = get_k2targets_table(select=select, cache=False)
        nasa_table = clean_kepler_table(nasa_table)
        
        # join info about planet host status: 
        full_k2cand_table = get_k2candidates_table()
        k2cand_table = unique(full_k2cand_table, keys='epic_name')
        k2cand_table_to_join = k2cand_table['k2c_disp','k2c_note']
        epic_numbers = [int(n.split(" ")[1]) for n in k2cand_table['epic_name']]
        k2cand_table_to_join['epic_number'] = epic_numbers
        nasa_table = join(nasa_table, k2cand_table_to_join, keys='epic_number', join_type='left')
                
        nasa_table_key = 'epic_number' # for join to nasa_table
        ra_key, dec_key = 'ra_epic', 'dec_epic' # for coordinate transforms
        ang_dist_key = 'k2_gaia_ang_dist' # name for gaia - nasa angular distance        
        
    elif exoplanets:
        gaia_matches_file = data_dir+'exoplanets_5arcsec_gaia.fits'
        dist_table_file = data_dir+'exoplanets_5arcsec_dist.fits'
        outfile_1arcsec = 'exoplanets_dr2_1arcsec.fits'
        
        select = None
        nasa_table = get_confirmed_planets_table(select=select, cache=False)
        nasa_table = clean_kepler_table(nasa_table)
        nasa_table.remove_columns(['ra','dec','ra_str','dec_str'])
        
                
        nasa_table_key = 'pl_name' # for join to nasa_table
        ra_key, dec_key = 'ra_nasa', 'dec_nasa' # for coordinate transforms
        ang_dist_key = 'nasa_gaia_ang_dist' # name for gaia - nasa angular distance
    else:
        print("one of the following flags must be True: kepler, k2, exoplanets")
        return
    
    hdus = fits.open(gaia_matches_file)
    gaia_matches_tbl = Table(hdus[1].data)
    gaia_matches_tbl = clean_gaia_table(gaia_matches_tbl, kepler=kepler, k2=k2, exoplanets=exoplanets)
    hdus = fits.open(dist_table_file)
    dist_tbl = Table(hdus[1].data)
    dist_tbl = clean_dist_table(dist_tbl)
    
    gaia_w_dist_tbl = join(gaia_matches_tbl, dist_tbl, keys='source_id', join_type='left')
    table = join(gaia_w_dist_tbl, nasa_table, keys=nasa_table_key)
        
    # calculate angular distances, propagating PM between epochs
    refCoord = coord.SkyCoord(ra=table[ra_key], dec=table[dec_key], obstime='J2000') # this is a guess!        
    table['radial_velocity'][np.isnan(table['radial_velocity'])] = 0.
    gaia_time = Time(table['gaia_ref_epoch'], format='jyear')
    ref_time = refCoord.obstime
    gaiaCoord = coord.SkyCoord(ra=table['ra'], 
                            dec=table['dec'], 
                            distance=(table['parallax']).to(u.pc, u.parallax()),
                            radial_velocity=table['radial_velocity'],
                            pm_ra_cosdec=table['pmra'], 
                            pm_dec=table['pmdec'], 
                            obstime=gaia_time
                            )
    gaiaCoord_shifted = gaiaCoord.apply_space_motion(new_obstime=ref_time)
    sep = refCoord.separation(gaiaCoord_shifted)
    ind = np.where(sep > 10. * u.deg)[0]
    for i in ind:
        sep[i] = 180.*u.deg - sep[i] # HACK
    table[ang_dist_key] = sep.arcsec
    table[ang_dist_key].unit = u.arcsec  
    
    if k2:
        # save 20 arcsec radius:
        table = table[table[ang_dist_key] <= 20.]
        table.write(outfile_20arcsec, format='fits', overwrite=True)
        print('{0} stars with matches within 20 arcsec'.format(len(np.unique(table[nasa_table_key]))))        
    
    if not exoplanets:
        # cut down to 4 arcsec and save:
        table = table[table[ang_dist_key] <= 4.]
        table.write(outfile_4arcsec, format='fits', overwrite=True)
        print('{0} stars with matches within 4 arcsec'.format(len(np.unique(table[nasa_table_key]))))

    # cut down to 1 arcsec and save:
    table = table[table[ang_dist_key] <= 1.]
    table.write(outfile_1arcsec, format='fits', overwrite=True)
    print('{0} stars with matches within 1 arcsec'.format(len(np.unique(table[nasa_table_key]))))  
    
    
    

if __name__ == "__main__":

    # Kepler:
    if False:
        make_full_tables(kepler=True)
        print('Kepler finished')
        
    
    # K2:
    if True:
        make_full_tables(k2=True)
        print('K2 finished')
        
    
    # confirmed planets:
    if False:
        make_full_tables(exoplanets=True)
        print('exoplanets finished')    
    