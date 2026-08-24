from astropy.table import Table
from scipy.optimize import curve_fit

#from astropy.coordinates import Angle, SkyCoord, search_around_sky

import astropy.coordinates as coord
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

catalog_dir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Catalogs/'
outdir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Outputs/Seperation_Matching_Outputs/'

def calc_bkg_ratio(x,y,r_tol, a):
    '''
        calculate bkg and src counts, and its ratio for a certain r_tol
        a is in units of /deg2
    '''
    #r_tol = 0.2
    count_src = sum( y[x<r_tol] )
    count_bkg = 0.5 * a *r_tol**2.
    return r_tol,count_src,count_bkg, count_bkg/count_src

Swift_157M_catalog = Table.read(catalog_dir + 'Swift_BAT_157m_Catalog_FITS.fits')
Fermi_4FGL_DR4 = Table.read(catalog_dir + '4FGLDR4.fit', hdu=1)

Fermi_Point_Source_Mask = [x == ''.ljust(18) for x in Fermi_4FGL_DR4['Extended_Source_Name'].tolist()]
Fermi_Point_Sources = Fermi_4FGL_DR4[Fermi_Point_Source_Mask]

BAT_RA = Swift_157M_catalog['RA'].tolist()
BAT_DEC = Swift_157M_catalog['DEC'].tolist()
LAT_RA = Fermi_Point_Sources['RAJ2000'].tolist()
LAT_DEC = Fermi_Point_Sources['DEJ2000'].tolist()

bat_c = coord.SkyCoord(BAT_RA, BAT_DEC, unit=(u.deg, u.deg), frame='icrs')
fermi_c = coord.SkyCoord(LAT_RA, LAT_DEC, unit=(u.deg, u.deg), frame='icrs')

dr = 0.02*u.deg ; r_max = 2.0
lis_r = np.arange(0, r_max, dr.value ) *u.deg 
r_4fit = 0.4 # r range for fitting

lis_N = np.zeros( len(lis_r) )
lis_Nsub = []
for r in np.arange(0, r_max+dr.value, dr.value )*u.deg :
    idxc, idxcatalog, d2d, d3d = bat_c.search_around_sky(fermi_c, r)
    lis_Nsub += [len(idxc)]
lis_N = np.array([lis_Nsub[i+1]-lis_Nsub[i] for i in range(len(lis_r))])

d = np.append( lis_r.value, lis_N ) ; d =d.reshape(2, len(lis_r) ).T

plt.figure()
plt.show()
plt.xlabel('radius (deg)')
plt.ylabel('N')
plt.axvline(0.08, color='black', lw=1., ls='--',label=r'$r_{sep} = 0.08^{o}$')

#fig, axes = plt.subplots(1, 2)
#ax1, ax2 = axes
#plt.plot(lis_r+dr/2., lis_N,'o', label='total: dr={0}'.format(dr) )
start = 0 * lis_r.unit               # <Quantity 0. deg>
x = np.append(start, lis_r + dr)     # now all in deg
plt.step(x.value, np.append(0, lis_N))
#plt.step(np.append(0,lis_r+dr), np.append(0,lis_N), label='total: dr={0}'.format(dr) )

r_thresh = 0.08
bin_centers = lis_r.value + dr.value/2.0  
mask = bin_centers < r_thresh
total_sources = np.sum(lis_N[mask])

print(f"Total sources with separation < {r_thresh}°: {total_sources}")

def model_distance_profile(x_data,m):
    return m*x_data 

x,y = (lis_r).value, lis_N
popt, pcov = curve_fit(model_distance_profile, x[x>0.2]*dr.value, y[x>0.2], p0=2.e3, sigma=None, absolute_sigma=False)
print(' a = {0:.3g} +/- {1:.3g} (/deg2)'.format( popt[0], pcov[0][0] ))
plt.plot(x, model_distance_profile(x*dr.value, popt[0]), ls='--', color='red', label=r'a={0:.2g}'.format( popt[0], pcov[0] ) )
plt.xlabel('(deg)')
plt.ylabel('#4FGL sources')
plt.legend()
plt.savefig(outdir + "Search_around_sky_result_hist.png",bbox_inches='tight')
plt.show()
plt.close()

lis_bkg_ratio = np.array( [ calc_bkg_ratio(x,y, r_tol, popt[0] )  for r_tol in np.arange(0,0.4, dr.value)  ] )

plt.plot( lis_bkg_ratio[:,0], lis_bkg_ratio[:,1], 'o-',label='src count')
plt.plot( lis_bkg_ratio[:,0], lis_bkg_ratio[:,2], 'o-',label='bkg count')
plt.plot( lis_bkg_ratio[:,0], lis_bkg_ratio[:,3]*100., 'o-',label='bkg ratio (%)')
plt.axhline(y=5.,c='grey',lw=1.,label='5% contamination')
plt.axvline(x=0.065, color='black', lw=1., ls='--')
plt.axvline(x=0.08, color = 'red', lw=1, ls = '--')
#plt.ylim(0,10)
#plt.xlim(0,0.2)
plt.legend() 
plt.xlabel('radius (deg)')
plt.savefig(outdir + "Signal_to_noise_plt")
plt.show()

idx_fermi, idx_bat, sep2d, sep3d = bat_c.search_around_sky(fermi_c, r_thresh *u.deg)

Swift_Matched = Swift_157M_catalog[idx_bat]
Fermi_Matched = Fermi_Point_Sources[idx_fermi]
Match_Seperation = [float(x.value) for x in sep2d]

Swift_Matched_Names = [Swift_Matched['COUNTERPART_NAME'][x] if Swift_Matched['COUNTERPART_NAME'][x] != 'NA' and Swift_Matched['COUNTERPART_NAME'][x] != 'None' else Swift_Matched['BAT_NAME'].tolist()[x] for x in range(len(Swift_Matched))]
Fermi_Matched_Names = [Fermi_Matched['ASSOC1'].tolist()[x].strip() if Fermi_Matched['ASSOC1'].tolist()[x].strip() != '' else Fermi_Matched['Source_Name'].tolist()[x].strip() for x in range(len(Fermi_Matched))]


Swift_Matched_RA = Swift_Matched['RA'].tolist()
Swift_Matched_DEC = Swift_Matched['DEC'].tolist()
Fermi_Matched_RA = Fermi_Matched['RAJ2000'].tolist()
Fermi_Matched_DEC = Fermi_Matched['DEJ2000'].tolist()

Swift_Matched_Class = Swift_Matched['TYPE'].tolist()
Fermi_Matched_Class = [x.strip() for x in Fermi_Matched['CLASS1'].tolist()]
Fermi_Matched_Class = ['unk' if x == '' else x for x in Fermi_Matched_Class]

Seperation_Match_Table = Table([Swift_Matched_Names, Fermi_Matched_Names, Match_Seperation, Swift_Matched_RA, Swift_Matched_DEC, Fermi_Matched_RA, Fermi_Matched_DEC, Swift_Matched_Class, Fermi_Matched_Class], names=['Swift_Name', 'Fermi_Name', 'Sep', 'Swift_RA', 'Swift_DEC', 'Fermi_RA', 'Fermi_DEC', 'Swift_Class', 'Fermi_Class'])
Seperation_Match_Table.write(outdir + 'Seperation_Matched_Table.fits', overwrite=True)
Seperation_Match_Table.pprint_all()