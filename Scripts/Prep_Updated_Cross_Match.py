from astropy.table import Table

import numpy as np

catalog_dir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Catalogs/'
output_dir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Outputs/'

# Load in Catalogs

Swift_157M_Catalog = Table.read(catalog_dir + 'Swift_BAT_157m_Catalog_FITS.fits')
Fermi_4FGL_Catalog = Table.read(catalog_dir + '4FGLDR4.fit', hdu=1)

Flagged_Matches = Table.read(output_dir + 'Update_Match_Flagging.fits')

# Mask and sort imported catalogs

Matched_Swift_Names = Flagged_Matches['Swift_Name'].tolist()
Matched_Fermi_Names = Flagged_Matches['Fermi_Name'].tolist()

Swift_Counterpart_Names = Swift_157M_Catalog['COUNTERPART_NAME'].tolist()
Swift_Source_Names = Swift_157M_Catalog['BAT_NAME'].tolist()


name_to_idx_swift = {}
for i, (src, cp) in enumerate(zip(Swift_Source_Names, Swift_Counterpart_Names)):
    name_to_idx_swift[src] = i
    name_to_idx_swift[cp] = i

matched_row_indices_swift = [name_to_idx_swift[name] for name in Matched_Swift_Names]

Matched_Swift_Sources = Swift_157M_Catalog[matched_row_indices_swift]

Fermi_Counterpart_Names = [x.strip() for x in Fermi_4FGL_Catalog['ASSOC1'].tolist()]
Fermi_Source_Names = [x.strip() for x in Fermi_4FGL_Catalog['Source_Name'].tolist()]

name_to_idx_fermi = {}
for i, (src, cp) in enumerate(zip(Fermi_Source_Names, Fermi_Counterpart_Names)):
    name_to_idx_fermi[src] = i
    name_to_idx_fermi[cp] = i

matched_row_indices_fermi = [name_to_idx_fermi[name] for name in Matched_Fermi_Names]

Matched_Fermi_Sources = Fermi_4FGL_Catalog[matched_row_indices_fermi]

# Get data from Swift Catalog

Matched_Swift_Gamma = Matched_Swift_Sources['GAMMA'].tolist()
Matched_Swift_Gamma_Err = Matched_Swift_Sources['GAMMA_ERR'].tolist()
Matched_Swift_Gamma_Errm, Matched_Swift_Gamma_Errp = zip(*[list(map(float, val.split('-'))) for val in Matched_Swift_Gamma_Err])

Matched_Swift_Flux = [float(x)*1e-12 for x in Matched_Swift_Sources['FLUX'].tolist()]
Matched_Swift_Flux_Err = Matched_Swift_Sources['FLUX_ERR'].tolist()
Matched_Swift_Flux_Errm, Matched_Swift_Flux_Errp = zip(*[list(map(float, val.split('-'))) for val in Matched_Swift_Flux_Err])

Matched_Swift_Flux_Errm = [x*1e-12 for x in Matched_Swift_Flux_Errm]
Matched_Swift_Flux_Errp = [x*1e-12 for x in Matched_Swift_Flux_Errp]

Matched_Swift_Flux = [list(np.array(x)*1e-3) for x in Matched_Swift_Sources['Spectral_Flux'].tolist()]
Matched_Swift_Flux_Err = [list(np.array(x)*1e-3) for x in Matched_Swift_Sources['Spectral_Flux_err'].tolist()]

# Get data from Fermi Catalog

Matched_Fermi_Gamma = Matched_Fermi_Sources['PL_Index'].tolist()
Matched_Fermi_Gamma_Err = [x if x is not None else 0 for x  in Matched_Fermi_Sources['Unc_PL_Index'].tolist()]

Matched_Fermi_Flux = Matched_Fermi_Sources['Energy_Flux100'].tolist()
Matched_Fermi_Flux_Err = Matched_Fermi_Sources['Unc_Energy_Flux100'].tolist()

Matched_Fermi_Spectral_Data = [[i* 624151 for i in x] for x in Matched_Fermi_Sources['nuFnu_Band'].tolist()]

Matched_Fermi_Flux_Data = Matched_Fermi_Sources['Flux_Band'].tolist()

Matched_Fermi_Flux_Data_Err = [list(zip(*x)) for x in Matched_Fermi_Sources['Unc_Flux_Band'].tolist()]
Matched_Fermi_Flux_Data_errm = [[abs(i) if i is not None else 0 for i in x[0]] for x in Matched_Fermi_Flux_Data_Err]
Matched_Fermi_Flux_Data_errp = [[abs(i) if i is not None else 0 for i in x[1]] for x in Matched_Fermi_Flux_Data_Err]

Matched_Fermi_Flux_errm_ratio = [np.array(a)/np.array(b) for a,b in list(zip(Matched_Fermi_Flux_Data_errm, Matched_Fermi_Flux_Data))]
Matched_Fermi_Flux_errp_ratio = [np.array(a)/np.array(b) for a,b in list(zip(Matched_Fermi_Flux_Data_errp, Matched_Fermi_Flux_Data))]

Matched_Fermi_Spectral_Data_errm = [list(np.array(a)*np.array(b)) for a,b in list(zip(Matched_Fermi_Flux_errm_ratio, Matched_Fermi_Spectral_Data))]
Matched_Fermi_Spectral_Data_errp = [list(np.array(a)*np.array(b)) for a,b in list(zip(Matched_Fermi_Flux_errp_ratio, Matched_Fermi_Spectral_Data))]

# Assemble SED Data

Cross_Matched_SED_Data = [x + y for x,y in list(zip(Matched_Swift_Flux, Matched_Fermi_Spectral_Data))]

Cross_Matched_SED_Data_errm = [x + y for x,y in list(zip(Matched_Swift_Flux_Err, Matched_Fermi_Spectral_Data_errm))]
Cross_Matched_SED_Data_errp = [x + y for x,y in list(zip(Matched_Swift_Flux_Err, Matched_Fermi_Spectral_Data_errp))]
Cross_Matched_SED_Data_err = list(zip(Cross_Matched_SED_Data_errm, Cross_Matched_SED_Data_errp))

Energy_Bin_Centers = [[0.016999999999999998, 0.022, 0.0295, 0.042499999999999996, 0.0625, 0.0875, 0.125, 0.1725, 70.71067811865474, 173.2050807568877, 547.722557505166, 1732.050807568877, 5477.225575051661, 17320.508075688773, 54772.25575051661, 316227.7660168379] for x in range(len(Matched_Fermi_Sources))]

# Assemble Catalog

Flagged_Match_Columns = [Flagged_Matches[x].tolist() for x in Flagged_Matches.colnames]
ids = np.arange(1, len(Flagged_Matches)+1, 1)

Cross_Match_Update_Catalog = Table([ids, Flagged_Match_Columns[0], Flagged_Match_Columns[1], Flagged_Match_Columns[3], Flagged_Match_Columns[2], Flagged_Match_Columns[4], Flagged_Match_Columns[9], Matched_Swift_Gamma, Matched_Swift_Gamma_Errm, Matched_Swift_Gamma_Errp, Matched_Swift_Flux, Matched_Swift_Flux_Errm, Matched_Swift_Flux_Errp, Flagged_Match_Columns[5], Flagged_Match_Columns[6], Matched_Fermi_Gamma, Matched_Fermi_Gamma_Err, Matched_Fermi_Flux, Matched_Fermi_Flux_Err, Flagged_Match_Columns[7], Flagged_Match_Columns[8], Energy_Bin_Centers, Cross_Matched_SED_Data, Cross_Matched_SED_Data_err], names=['id', 'flag', 'bat_name_counterpart', 'bat_category_type', 'fermi_name_counterpart', 'fermi_category_type', 'separation', 'bat_pindex', 'bat_pindex_errm', 'bat_pindex_errp', 'bat_flux', 'bat_flux_errm', 'bat_flux_errp', 'bat_ra', 'bat_dec', 'fermi_pindex_PL', 'fermi_pindex_PL_err','fermi_flux', 'fermi_flux_err' ,'fermi_ra', 'fermi_dec', 'Energy', 'Energy_Flux', 'Energy_Flux_Err'])
Cross_Match_Update_Catalog.write(output_dir + 'cross_match_updated.fits', overwrite=True)
