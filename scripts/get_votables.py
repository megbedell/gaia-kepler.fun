import os
from astropy.utils.data import download_file

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
                      
def save_table(url, out_name, select=None, format='votable', show_progress=True):
    url += '&format='+format
    if not select is None:
        url += '&select='+select
    table_path = download_file(url, cache=False,
                                   show_progress=show_progress,
                                   timeout=300)
    os.rename(table_path, out_name)
    
def kepler_table(out_name, format='votable', show_progress=True):
    select = 'kepid,tm_designation,ra,dec,kepmag'
    select += ',teff,teff_err1,teff_err2,teff_prov,logg,logg_err1,logg_err2,logg_prov'
    select += ',feh,feh_err1,feh_err2,feh_prov,radius,radius_err1,radius_err2'
    select += ',mass,mass_err1,mass_err2,prov_sec,nconfp,nkoi,ntce,jmag,hmag,kmag'
    save_table(STELLAR_URL, out_name, select=select, format=format, 
                                   show_progress=show_progress)
    print('Kepler table downloaded and saved as:')
    print(out_name)
    
def k2_table(out_name, format='votable', show_progress=True):
    select = 'epic_number,tm_name,k2_campaign_str,k2_type,ra,dec,k2_lcflag,k2_scflag'
    select += ',k2_teff,k2_tefferr1,k2_tefferr2,k2_logg,k2_loggerr1,k2_loggerr2'
    select += ',k2_metfe,k2_metfeerr1,k2_metfeerr2,k2_rad,k2_raderr1,k2_raderr2'
    select += ',k2_mass,k2_masserr1,k2_masserr2'
    save_table(K2TARGETS_URL, out_name, select=select, format=format, 
                                   show_progress=show_progress)
    print('K2 table downloaded and saved as:')
    print(out_name)
    
def confirmed_table(out_name, format='votable', show_progress=True):
    save_table(EXOPLANETS_URL, out_name, select=None, format=format, 
                                   show_progress=show_progress)
    print('Confirmed planets table downloaded and saved as:')
    print(out_name)
    
                      
if __name__ == "__main__":
    kepler_table('../data/keplerstellar.votable')
    k2_table('../data/k2targets.votable')
    confirmed_table('../data/exoplanets.votable')