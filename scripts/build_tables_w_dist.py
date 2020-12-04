from astropy.table import QTable, Table, Column, MaskedColumn, join, unique, vstack
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits, ascii
import astropy.coordinates as coord
from astropy.time import Time
from units import gaia_unit_map, kepler_unit_map
from nasa_tables import *
from pyia import GaiaData 
import pdb

no_dist = True # flag to not include Bailer-Jones distances

def clean_gaia_table(tbl, kepler=False, k2=False,exoplanets=False):
    """
    Add units, delete some columns
    """
    '''
    cols_to_delete = ['solution_id', 'random_index', 'astrometric_n_obs_al', 'astrometric_n_obs_ac', 
                      'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al',
                      'astrometric_params_solved', 'astrometric_weight_al', 'astrometric_pseudo_colour',
                      'astrometric_pseudo_colour_error', 'mean_varpi_factor_al', 'astrometric_matched_observations',
                      'visibility_periods_used', 'astrometric_sigma5d_max',  
                      'frame_rotator_object_type', 'matched_observations', 'phot_g_n_obs', 
                      'phot_g_mean_flux_over_error', 'phot_bp_n_obs', 'phot_bp_mean_flux_over_error',
                      'phot_rp_n_obs', 'phot_rp_mean_flux_over_error', 
                      'phot_proc_mode', 'rv_nb_transits', 'rv_template_teff', 'rv_template_logg', 
                      'rv_template_fe_h', 'priam_flags', 'flame_flags', 'datalink_url', 'epoch_photometry_url',
                      'phot_bp_n_contaminated_transits', 'phot_bp_n_blended_transits',
                      'phot_rp_n_contaminated_transits', 'phot_rp_n_blended_transits',
                      'scan_direction_strength_k1', 'scan_direction_mean_k1',
                      'scan_direction_strength_k2', 'scan_direction_mean_k2',
                      'scan_direction_strength_k3', 'scan_direction_mean_k3',
                      'scan_direction_strength_k4', 'scan_direction_mean_k4',
                      'ipd_gof_harmonic_amplitude', 'ipd_gof_harmonic_phase', 
                      'ipd_frac_multi_peak', 'ipd_frac_odd_win'
                    ]
    tbl.remove_columns(cols_to_delete)
    '''
    if kepler:
        #cols_to_delete = ['kepler_oid', 'angdist']
        cols_to_delete = ['kepler_oid']
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
        gaia_matches_file = data_dir+'kepler_30arcsec_gaia.fits'
        dist_table_file = data_dir+'kepler_30arcsec_dist.fits'
        outfile_20arcsec = 'kepler_dr2_20arcsec.fits'
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
        gaia_matches_file = data_dir+'k2_30arcsec_gaia.fits'
        dist_table_file = data_dir+'k2_30arcsec_dist.fits'
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
        gaia_matches_file = data_dir+'exoplanets_10arcsec_gaia.fits'
        dist_table_file = data_dir+'exoplanets_10arcsec_dist.fits'
        outfile_1arcsec = 'exoplanets_dr2_1arcsec.fits'
        
        select = 'pl_hostname,pl_letter,pl_name,pl_discmethod,pl_controvflag,pl_pnum,pl_orbper'
        select += ',pl_orbpererr1,pl_orbpererr2,pl_orbperlim,pl_orbpern,pl_orbsmax,pl_orbsmaxerr1,pl_orbsmaxerr2'
        select += ',pl_orbsmaxlim,pl_orbsmaxn,pl_orbeccen,pl_orbeccenerr1,pl_orbeccenerr2,pl_orbeccenlim'
        select += ',pl_orbeccenn,pl_orbincl,pl_orbinclerr1,pl_orbinclerr2,pl_orbincllim,pl_orbincln'
        select += ',pl_bmassj,pl_bmassjerr1,pl_bmassjerr2,pl_bmassjlim,pl_bmassn,pl_bmassprov'
        select += ',pl_radj,pl_radjerr1,pl_radjerr2,pl_radjlim,pl_radn,pl_dens,pl_denserr1'
        select += ',pl_denserr2,pl_denslim,pl_densn,pl_ttvflag,pl_kepflag,pl_k2flag'
        select += ',st_posn,st_dist,st_disterr1'
        select += ',st_disterr2,st_distlim,st_distn,st_optmag,st_optmagerr,st_optmaglim,st_optband'
        select += ',gaia_gmag,gaia_gmagerr,gaia_gmaglim,st_teff,st_tefferr1,st_tefferr2,st_tefflim'
        select += ',st_teffn,st_mass,st_masserr1,st_masserr2,st_masslim,st_massn,st_rad,st_raderr1'
        select += ',st_raderr2,st_radlim,st_radn,pl_nnotes,rowupdate,pl_facility'
        select += ',pl_tranflag,pl_rvflag,pl_imgflag,pl_astflag,pl_omflag,pl_cbflag,pl_disc'
        select += ',hd_name,hip_name,st_radv,st_radverr1,st_radverr2,st_radvlim,st_radvn'
        
        nasa_table = get_confirmed_planets_table(select=select, cache=False)
        nasa_table = clean_kepler_table(nasa_table)
        #nasa_table.remove_columns(['ra','dec','ra_str','dec_str'])
        
                
        nasa_table_key = 'pl_name' # for join to nasa_table
        ra_key, dec_key = 'ra_nasa', 'dec_nasa' # for coordinate transforms
        ang_dist_key = 'nasa_gaia_ang_dist' # name for gaia - nasa angular distance
    else:
        print("one of the following flags must be True: kepler, k2, exoplanets")
        return
    
    hdus = fits.open(gaia_matches_file)
    gaia_matches_tbl = Table(hdus[1].data)
    gaia_matches_tbl = clean_gaia_table(gaia_matches_tbl, kepler=kepler, k2=k2, exoplanets=exoplanets)
    if no_dist:
        table = join(gaia_matches_tbl, nasa_table, keys=nasa_table_key)
    else:
        hdus = fits.open(dist_table_file)
        dist_tbl = Table(hdus[1].data)
        dist_tbl = clean_dist_table(dist_tbl)
        gaia_w_dist_tbl = join(gaia_matches_tbl, dist_tbl, keys='source_id', join_type='left')
        table = join(gaia_w_dist_tbl, nasa_table, keys=nasa_table_key)
        
    # calculate angular distances, propagating PM between epochs
    table.sort('source_id')
    coordtable = GaiaData(gaia_matches_file)
    coordtable.data.sort('source_id')
    gaiaCoord = coordtable.get_skycoord(distance=coordtable.get_distance(min_parallax=1e-3*u.mas,
                                    parallax_fill_value=1e-5*u.mas))
    # NOTE: the above does not use RVs at all. will need to deal with EDR3 naming issue.
       
    refCoord = coord.SkyCoord(ra=table[ra_key], dec=table[dec_key], 
                            obstime='J2000') # KIC/EPIC/NExScI coord epoch
    ref_time = refCoord.obstime

    gaiaCoord_shifted = gaiaCoord.apply_space_motion(new_obstime=ref_time)
    sep = refCoord.separation(gaiaCoord_shifted)
    pdb.set_trace()
    ind = np.where(sep > 10. * u.deg)[0]
    for i in ind:
        sep[i] = 180.*u.deg - sep[i] # HACK
    table[ang_dist_key] = sep.arcsec
    table[ang_dist_key].unit = u.arcsec 
    table.sort([nasa_table_key, ang_dist_key]) # within each KIC/EPIC/host, sort by ang sep

    if not exoplanets:
        # save 20 arcsec radius:
        table = table[table[ang_dist_key] <= 20.]
        table.write(outfile_20arcsec, format='fits', overwrite=True)
        print('{0} stars with matches within 20 arcsec'.format(len(np.unique(table[nasa_table_key]))))
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
    if True:
        make_full_tables(kepler=True)
        print('Kepler finished')
        
    
    # K2:
    if False:
        make_full_tables(k2=True)
        print('K2 finished')
        
    
    # confirmed planets:
    if False:
        make_full_tables(exoplanets=True)
        print('exoplanets finished')    
    