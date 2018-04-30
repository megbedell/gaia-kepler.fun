from astropy.table import QTable, Table, Column, MaskedColumn, join, unique, vstack
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
import astropy.coordinates as coord
from astropy.time import Time
from units import gaia_unit_map, kepler_unit_map
from nasa_tables import *

def clean_gaia_table(tbl, kepler=False, k2=False):
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
        cols_to_delete = ['kepler_oid', 'ra_kic', 'dec_kic', 'angdist', 'kepid']
        tbl.remove_columns(cols_to_delete)
    if k2:
        cols_to_delete = ['k2_oid', 'ra_epic', 'dec_epic', 'dist', 'epic_number']
        tbl.remove_columns(cols_to_delete)
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
                
    tbl.rename_column('ra', 'ra_kepler')
    tbl.rename_column('dec', 'dec_kepler')
    return tbl
    

if __name__ == "__main__":

    # KEPLER:
    if True:
        data_dir = '../data/'
        kepler_gaia_file = data_dir+'kepler_4arcsec.fits'
        hdus = fits.open(kepler_gaia_file)
        kepler_gaia_tbl = Table(hdus[1].data)
        kepler_gaia_tbl = clean_gaia_table(kepler_gaia_tbl, kepler=True)
    
        kepler_dist_file = data_dir+'kepler_4arcsec_w_dist.fits'
        hdus = fits.open(kepler_dist_file)
        kepler_dist_tbl = Table(hdus[1].data)
        kepler_dist_tbl = clean_dist_table(kepler_dist_tbl)

        select = 'kepid,tm_designation,ra,dec,kepmag'
        select += ',teff,teff_err1,teff_err2,teff_prov,logg,logg_err1,logg_err2,logg_prov'
        select += ',feh,feh_err1,feh_err2,feh_prov,radius,radius_err1,radius_err2'
        select += ',mass,mass_err1,mass_err2,prov_sec,nconfp,nkoi,ntce,jmag,hmag,kmag'
        kepler_nasa_table = get_keplerstellar_table(select=select)
        kepler_nasa_table = clean_kepler_table(kepler_nasa_table)
    
        kepler_gaia_w_dist_tbl = join(kepler_gaia_tbl, kepler_dist_tbl, keys='source_id')
        table = join(kepler_gaia_w_dist_tbl, kepler_nasa_table, keys='kepid')
        table['planet?'] = 'none'
        table['planet?'][table['nkoi'] > 0] = 'cand'
        table['planet?'][table['nconfp'] > 0] = 'conf'
    
        kicCoord = coord.SkyCoord(ra=table['ra_kepler'], dec=table['dec_kepler'])
        table['radial_velocity'][np.isnan(table['radial_velocity'])] = 0.
        gaia_time = Time(table['gaia_ref_epoch'], format='jyear')
        gaiaCoord = coord.SkyCoord(ra=table['ra'], 
                                dec=table['dec'], 
                                distance=(table['parallax']).to(u.pc, u.parallax()),
                                radial_velocity=table['radial_velocity'],
                                pm_ra_cosdec=table['pmra'], 
                                pm_dec=table['pmdec'], 
                                obstime=gaia_time
                                )
        sep = kicCoord.separation(gaiaCoord)
        ind = np.where(sep > 10. * u.deg)[0]
        for i in ind:
            sep[i] = 180.*u.deg - sep[i] # HACK
        table['kepler_gaia_ang_dist'] = sep.arcsec
        table['kepler_gaia_ang_dist'].unit = u.arcsec
        table.write('kepler_dr2_4arcsec.fits', format='fits', overwrite=True)
        print('{0} stars with matches within 4 arcsec'.format(len(np.unique(table['kepid']))))
    
    
        table = table[table['kepler_gaia_ang_dist'] <= 1.]
        table.write('kepler_dr2_1arcsec.fits', format='fits', overwrite=True)
    
        print('{0} stars with matches within 1 arcsec'.format(len(np.unique(table['kepid']))))
        print('Kepler finished')
        
    
    # K2:
    if True:
        data_dir = '../data/'
        k2_gaia_file = data_dir+'k2_4arcsec.fits'
        hdus = fits.open(k2_gaia_file)
        k2_gaia_tbl = Table(hdus[1].data)
        k2_gaia_tbl = clean_gaia_table(k2_gaia_tbl, k2=True)
    
        k2_dist_file = data_dir+'k2_4arcsec_w_dist.fits'
        hdus = fits.open(k2_dist_file)
        k2_dist_tbl = Table(hdus[1].data)
        k2_dist_tbl = clean_dist_table(k2_dist_tbl)
    
        select = 'epic_number,tm_name,k2_campaign_str,k2_type,ra,dec,k2_lcflag,k2_scflag'
        select += ',k2_teff,k2_tefferr1,k2_tefferr2,k2_logg,k2_loggerr1,k2_loggerr2'
        select += ',k2_metfe,k2_metfeerr1,k2_metfeerr2,k2_rad,k2_raderr1,k2_raderr2'
        select += ',k2_mass,k2_masserr1,k2_masserr2'
        k2_nasa_table = get_k2targets_table(select=select)
        k2_nasa_table = clean_kepler_table(k2_nasa_table)
    
        k2_gaia_w_dist_tbl = join(k2_gaia_tbl, k2_dist_tbl, keys='source_id')
        table = join(k2_gaia_w_dist_tbl, k2_nasa_table, keys='epic_number')
    
        full_k2cand_table = get_k2candidates_table()
        k2cand_table = unique(full_k2cand_table, keys='epic_name')
        k2cand_table_to_join = k2cand_table['k2c_disp','k2c_note']
        epic_numbers = [int(n.split(" ")[1]) for n in k2cand_table['epic_name']]
        k2cand_table_to_join['epic_number'] = epic_numbers
        table = join(table, k2cand_table_to_join, keys='epic_number', join_type='left')
    
        epicCoord = coord.SkyCoord(ra=table['ra_kepler'], dec=table['dec_kepler'])
        table['radial_velocity'][np.isnan(table['radial_velocity'])] = 0.
        gaia_time = Time(table['gaia_ref_epoch'], format='jyear')
        gaiaCoord = coord.SkyCoord(ra=table['ra'], 
                                dec=table['dec'], 
                                distance=(table['parallax']).to(u.pc, u.parallax()),
                                radial_velocity=table['radial_velocity'],
                                pm_ra_cosdec=table['pmra'], 
                                pm_dec=table['pmdec'], 
                                obstime=gaia_time
                                )
        sep = epicCoord.separation(gaiaCoord)
        ind = np.where(sep > 10. * u.deg)[0]
        for i in ind:
            sep[i] = 180.*u.deg - sep[i] # HACK
        table['kepler_gaia_ang_dist'] = sep.arcsec
        table['kepler_gaia_ang_dist'].unit = u.arcsec
        table.write('k2_dr2_4arcsec.fits', format='fits', overwrite=True)
        print('{0} stars with matches within 4 arcsec'.format(len(np.unique(table['epic_number']))))
        
    
        table = table[table['kepler_gaia_ang_dist'] <= 1.]
        table.write('k2_dr2_1arcsec.fits', format='fits', overwrite=True)
    
        print('{0} stars with matches within 1 arcsec'.format(len(np.unique(table['epic_number']))))
        print('K2 finished')
        
    
    # confirmed planets:
    #conf_nasa_table = get_confirmed_planets_table(select=None)
    
    