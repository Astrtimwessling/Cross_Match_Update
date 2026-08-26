from astroquery.simbad import Simbad
from astropy.table import Table
from collections import defaultdict
from astropy.coordinates import Angle

import astropy.coordinates as coord
import astropy.units as u
import numpy as np

catalog_dir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Catalogs/'
outdir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Outputs/Identification_Matching_Outputs/'

def resolve_to_main_id(name_list):
    """
    Query SIMBAD for each name and return a dict mapping
    the original name to the SIMBAD MAIN_ID.
    """
    result = custom_simbad.query_objects(name_list) 
    main_ids = result['main_id']                   
    return dict(zip(name_list, main_ids))

Swift_157M_catalog = Table.read(catalog_dir + 'Swift_BAT_157m_Catalog_FITS.fits')
Fermi_4FGL_DR4 = Table.read(catalog_dir + '4FGLDR4.fit', hdu=1)

Fermi_Point_Source_Mask = [x == ''.ljust(18) for x in Fermi_4FGL_DR4['Extended_Source_Name'].tolist()]
Fermi_Point_Sources = Fermi_4FGL_DR4[Fermi_Point_Source_Mask]

custom_simbad = Simbad()
custom_simbad.reset_votable_fields()             
custom_simbad.add_votable_fields('ids')

Fermi_Names = Fermi_Point_Sources['ASSOC1'].tolist()
Swift_Names = Swift_157M_catalog['COUNTERPART_NAME'].tolist()

Fermi_Names = [x.replace('[HB89]', '[HB93]') if "[HB89]" in x else x.replace("MG4", "MITG") if "MG4" in x else x for x in Fermi_Names]
Swift_Names = [x.replace('[HB89]', '[HB93]') if "[HB89]" in x else x.replace("MG4", "MITG") if "MG4" in x else x for x in Swift_Names]

fermi_map = resolve_to_main_id(Fermi_Names)
swift_map = resolve_to_main_id(Swift_Names)

fermi_set = {mid for mid in fermi_map.values() if mid is not None}
swift_set = {mid for mid in swift_map.values() if mid is not None}  

common_ids = swift_set & fermi_set 

rev_swift = defaultdict(list)
for orig, mid in swift_map.items():
    if mid in common_ids:
        rev_swift[mid].append(orig)

rev_fermi = defaultdict(list)
for orig, mid in fermi_map.items():
    if mid in common_ids:
        rev_fermi[mid].append(orig)
        
simbad_swift_matched_names = []
simbad_fermi_matched_names = []

for mid in common_ids:
    simbad_swift_matched_names.append(rev_swift[mid])
    simbad_fermi_matched_names.append(rev_fermi[mid])
    
simbad_swift_matched_names.remove(simbad_swift_matched_names[0])
simbad_fermi_matched_names.remove(simbad_fermi_matched_names[0])

simbad_swift_matched_names = [x[0] for x in simbad_swift_matched_names]
simbad_fermi_matched_names = [x[0].strip() for x in simbad_fermi_matched_names]

id_matched_swift_names = [x.replace('[HB93]', '[HB89]',) if '[HB93]' in x else x.replace("MITG", "MG4") if 'MITG' in x else x for x in simbad_swift_matched_names] + ['SGR A*', "Kepler's SNR", 'B2 2023+33', 'B2 0920+39']
id_matched_fermi_names = [x.replace('[HB93]', '[HB89]',) if '[HB93]' in x else x.replace("MITG", "MG4") if 'MITG' in x else x for x in simbad_fermi_matched_names] + ['Galactic Center', "Kepler SNR", 'B2 2023+33', 'B2 0920+39']

print(f"Total id matches found: {len(id_matched_swift_names)}")

idx_order_swift = [{name: i for i, name in enumerate(Swift_157M_catalog['COUNTERPART_NAME'].tolist())}[name] for name in id_matched_swift_names]
Swift_Matched = Swift_157M_catalog[idx_order_swift]

idx_order_fermi = [{name: i for i, name in enumerate([x.strip() for x in Fermi_Point_Sources['ASSOC1'].tolist()])}[name] for name in id_matched_fermi_names]
Fermi_Matched = Fermi_Point_Sources[idx_order_fermi]

Swift_Matched_Names = [Swift_Matched['COUNTERPART_NAME'][x] if Swift_Matched['COUNTERPART_NAME'][x] != 'NA' and Swift_Matched['COUNTERPART_NAME'][x] != 'None' else Swift_Matched['BAT_NAME'].tolist()[x] for x in range(len(Swift_Matched))]
Fermi_Matched_Names = [Fermi_Matched['ASSOC1'].tolist()[x].strip() if Fermi_Matched['ASSOC1'].tolist()[x].strip() != '' else Fermi_Matched['ASSOC2'].tolist()[x].strip() if Fermi_Matched['ASSOC2'].tolist()[x].strip() != '' else Fermi_Matched['Source_Name'].tolist()[x].strip() for x in range(len(Fermi_Matched))]

Swift_Matched_RA = Swift_Matched['RA'].tolist()
Swift_Matched_DEC = Swift_Matched['DEC'].tolist()
Fermi_Matched_RA = Fermi_Matched['RAJ2000'].tolist()
Fermi_Matched_DEC = Fermi_Matched['DEJ2000'].tolist()

swift_c = coord.SkyCoord(np.array(Swift_Matched_RA), np.array(Swift_Matched_DEC), unit=(u.deg, u.deg), frame='icrs')
fermi_c = coord.SkyCoord(np.array(Fermi_Matched_RA), np.array(Fermi_Matched_DEC), unit=(u.deg, u.deg), frame='icrs')

Match_Seperation = list(Angle(swift_c.separation(fermi_c), unit=u.degree).degree)

Swift_Matched_Class = Swift_Matched['TYPE'].tolist()
Fermi_Matched_Class = Fermi_Matched['CLASS1'].tolist()
Fermi_Matched_Class = ['unk' if x == '' else x for x in Fermi_Matched_Class]

Identification_Match_Table = Table([Swift_Matched_Names, Fermi_Matched_Names, Match_Seperation, Swift_Matched_RA, Swift_Matched_DEC, Fermi_Matched_RA, Fermi_Matched_DEC, Swift_Matched_Class, Fermi_Matched_Class], names=['Swift_Name', 'Fermi_Name', 'Sep', 'Swift_RA', 'Swift_DEC', 'Fermi_RA', 'Fermi_DEC', 'Swift_Class', 'Fermi_Class'])
Identification_Match_Table.write(outdir+'Identification_Match_Table.fits',overwrite=True)