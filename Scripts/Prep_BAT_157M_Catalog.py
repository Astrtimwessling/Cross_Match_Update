from astropy.table import Table, hstack

import numpy as np

catalog_dir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Catalogs/'

SED_Data = Table.read(catalog_dir + 'BAT_157m_catalog_SED.fits')
Swift_157M_catalog = Table.read(catalog_dir + "Swift_BAT_157m.csv", format="csv")

Source_Names = Swift_157M_catalog['BAT Name (a)'].tolist()
Source_Names = [s.replace("SWIFT ", "SWIFT_") for s in Source_Names]
Source_Names = [s.replace("Swift ", "SWIFT_") for s in Source_Names]

Source_id = Swift_157M_catalog['#'].tolist()

SED_data_names = SED_Data['BAT_name'].tolist()
SED_data_names = [s.replace("Swift ", "SWIFT_") for s in SED_data_names]

order_map = {val: i for i, val in enumerate(Source_Names)}
sort_idx = np.argsort([order_map[val] for val in SED_data_names])

SED_data_ordered = SED_Data[sort_idx]
SED_data_ordered.remove_column('BAT_name')

Source_Names_Table = Table([Source_id, Source_Names], names=['Index', 'BAT_name'])

Swift_157M_catalog_FITS = hstack([Source_Names_Table,Swift_157M_catalog,SED_data_ordered])
Swift_157M_catalog_FITS.remove_columns(['BAT Name (a)', '#'])

Swift_157M_catalog_FITS.rename_columns(Swift_157M_catalog_FITS.colnames ,['Index', 'BAT_NAME', 'RA', 'DEC', 'SNR', 'COUNTERPART_NAME', 'CTPT_RA', 'CTPT_DEC', 'FLUX', 'FLUX_ERR', 'GAMMA', 'GAMMA_ERR', 'CHI_SQ_R', 'REDSHIFT', 'LUM', 'ASSOC_STERN', 'CL2', 'TYPE', 'Energy', 'Energy_err', 'Spectral_Flux', 'Spectral_Flux_err'])

Swift_157M_catalog_FITS['RA'].unit = 'deg'
Swift_157M_catalog_FITS['DEC'].unit = 'deg'
Swift_157M_catalog_FITS['CTPT_RA'].unit = 'deg'
Swift_157M_catalog_FITS['CTPT_DEC'].unit = 'deg'
Swift_157M_catalog_FITS['FLUX'].unit = '10^-12er'
Swift_157M_catalog_FITS['FLUX_ERR'].unit = '10^-12er'
Swift_157M_catalog_FITS['Energy'].unit = 'keV'
Swift_157M_catalog_FITS['Energy_err'].unit = 'keV'
Swift_157M_catalog_FITS['Spectral_Flux'].unit = 'keV / (s cm2)'
Swift_157M_catalog_FITS['Spectral_Flux_err'].unit = 'keV / (s cm2)'

Swift_157M_catalog_FITS.write(catalog_dir + 'Swift_BAT_157m_Catalog_FITS.fits', overwrite=True)




