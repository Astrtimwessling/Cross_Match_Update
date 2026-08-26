from astropy.coordinates import Angle, SkyCoord

import astropy.coordinates as coord
import astropy.units as u

def source_separation(fermi_name, swift_name, fermi_catalog, swift_catalog):
    fermi_ASSOC1_names = [x.strip() for x in fermi_catalog['ASSOC1'].tolist()]
    fermi_ASSOC2_names = [x.strip() for x in fermi_catalog['ASSOC2'].tolist()]
    fermi_source_names = [x.strip() for x in fermi_catalog['Source_Name'].tolist()]
    
    if fermi_name in fermi_ASSOC1_names:
        idx = fermi_ASSOC1_names.index(fermi_name)
        
    elif fermi_name in fermi_ASSOC2_names:
        idx = fermi_ASSOC2_names.index(fermi_name)
        
    elif fermi_name in fermi_source_names:
        idx = fermi_source_names.index(fermi_name)
        
    else:
        print(f"{fermi_name} not found in Fermi Catalog")
        
    source_entry_fermi = fermi_catalog[idx]
    
    fermi_RA = source_entry_fermi['RAJ2000'].tolist()
    fermi_DEC = source_entry_fermi['DEJ2000'].tolist()

    fermi_c = coord.SkyCoord(fermi_RA, fermi_DEC, unit=(u.deg, u.deg), frame='icrs')
    
    swift_counterpart_names = swift_catalog['COUNTERPART_NAME'].tolist()
    swift_source_names = swift_catalog['BAT_NAME'].tolist()
    
    if swift_name in swift_counterpart_names:
        idx2 = swift_counterpart_names.index(swift_name)
    elif swift_name in swift_source_names:
        idx2 = swift_source_names.index(swift_name)
        
    source_entry_swift = swift_catalog[idx2]

    swift_RA = source_entry_swift['RA'].tolist()
    swift_DEC = source_entry_swift['DEC'].tolist()
    
    swift_c = coord.SkyCoord(swift_RA, swift_DEC, unit=(u.deg, u.deg), frame='icrs')
    sep2d = swift_c.separation(fermi_c)

    angle = Angle(sep2d, unit=u.degree)
    sep2d_deg = angle.degree

    return sep2d_deg