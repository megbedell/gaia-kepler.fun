import numpy as np
import altair as alt
import pandas as pd
from astropy.table import Table

def prepare_data(table, subsample_size=4999):
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
    return subsample

def skyview_cmd(subsample, out_file, data_type='kepler'):
    brush = alt.selection(type='interval', resolve='global', 
                          on="[mousedown[event.shiftKey], window:mouseup] > \
                          window:mousemove!", zoom='False',
                          translate="[mousedown[event.shiftKey], window:mouseup] > \
                          window:mousemove!")

    pan = alt.selection(type='interval', bind='scales',
                        on="[mousedown[!event.shiftKey], window:mouseup] > \
                        window:mousemove!",
                        translate="[mousedown[!event.shiftKey], window:mouseup] > \
                        window:mousemove!")

    if data_type=='kepler':
        scale = alt.Scale(domain=['none', 'cand', 'conf'],
                      range=['#4a69bd', '#e55039', '#f6b93b'])
        color = alt.Color('planet?:N', scale=scale)
        xscale = alt.Scale(domain=[270, 310])
        yscale = alt.Scale(domain=[35, 55])
    elif data_type=='k2':
        #scale = alt.Scale(domain=['none', 'cand'],
        #              range=['#4a69bd', '#e55039'])
        #color = alt.Color('k2c_disp:N', scale=scale)
        color = 'k2c_disp'
        xscale = alt.Scale(domain=[0, 360])
        yscale = alt.Scale(domain=[-90, 90])
    else:
        print("ERROR: data_type not recognized")
        return

    chart1 = alt.Chart(subsample).mark_point().encode(
            alt.X('ra_gaia', scale=xscale, 
                  axis=alt.Axis(title='Gaia RA (deg)')),
            alt.Y('dec_gaia', scale=yscale,
                  axis=alt.Axis(title='Gaia Dec (deg)')),
            color=alt.condition(brush, color, 
                                alt.ColorValue('gray'))
        ).properties(
    selection=brush+pan,
    projection={'type': 'gnomonic'},
    width=450,
    height=500
    )

    chart2 = alt.Chart(subsample).mark_point().encode(
            alt.X('gaiamag_minus_kepmag', axis=alt.Axis(title='G - K (mag)')),
            alt.Y('abs_gmag', axis=alt.Axis(title='abs. G')),
            color=alt.condition(brush, color, 
                            alt.ColorValue('gray'))
        ).properties(
    selection=brush,
    width=450,
    height=500
    )

    chart = chart1 | chart2 
    chart.savechart(out_file)
    


if __name__ == "__main__":
    subsample_size = 4999
    
    kep_data = Table.read('../data/kepler_tgas_1arcsec.fits', format='fits')
    subsample = prepare_data(kep_data, subsample_size=subsample_size)
    skyview_cmd(subsample, '../kepler_tgas_demo.html', data_type='kepler')
    
    k2_data = Table.read('../data/k2_tgas_1arcsec.fits', format='fits')
    subsample = prepare_data(k2_data, subsample_size=subsample_size)
    skyview_cmd(subsample, '../k2_tgas_demo.html', data_type='k2')
 
