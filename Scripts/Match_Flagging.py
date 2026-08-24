from astropy.table import Table
from collections import Counter

import numpy as np

catalog_dir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Catalogs/'
output_dir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Outputs'

Spatial_Matches = Table.read(output_dir + '/Seperation_Matching_Outputs/Seperation_Matched_Table.fits')
Identification_Matches = Table.read(output_dir + '/Identification_Matching_Outputs/Identification_Match_Table.fits')

Minami_Flags = Table.read(catalog_dir + 'crossmatch_update+LAT_DR4model_260822.fits')

Spatial_Match_Names = list(zip(Spatial_Matches['Swift_Name'].tolist(), Spatial_Matches['Fermi_Name'].tolist()))
Identification_Match_Names = list(zip(Identification_Matches['Swift_Name'].tolist(), Identification_Matches['Fermi_Name'].tolist()))

Both_Matched = [x for x in Spatial_Match_Names if x in Identification_Match_Names]
Only_id_Matched = [x for x in Identification_Match_Names if x not in Spatial_Match_Names]
Only_spatial_Matched = [x for x in Spatial_Match_Names if x not in Identification_Match_Names]
All_Matches = Spatial_Match_Names + Only_id_Matched

print(f"Number spatial and id matched: {len(Both_Matched)}, Number only spatial matched: {len(Only_spatial_Matched)}, Number only id matched: {len(Only_id_Matched)}")
print(f"Total Matched: {len(Both_Matched) + len(Only_id_Matched) + len(Only_spatial_Matched)}")

Only_id_Matched_mask = [x in Identification_Match_Names and x not in Spatial_Match_Names for x in Identification_Match_Names]
Id_Only_Matches = Identification_Matches[Only_id_Matched_mask]
Id_Only_Matches = Id_Only_Matches[Id_Only_Matches.argsort('Swift_RA')]


fermi_types = list(Counter(Spatial_Matches['Fermi_Class'].tolist() + Id_Only_Matches['Fermi_Class'].tolist()).keys())
swift_types = list(Counter(Spatial_Matches['Swift_Class'].tolist() + Id_Only_Matches['Swift_Class'].tolist()).keys())

fermi_source_types = [x.strip() for x in Spatial_Matches['Fermi_Class'].tolist() + Id_Only_Matches['Fermi_Class'].tolist()]
swift_source_types = Spatial_Matches['Swift_Class'].tolist() + Id_Only_Matches['Swift_Class'].tolist()

Fermi_Blazars = ['BCU','BLL','FSRQ','fsrq','bll','bcu']
Swift_Blazars = ['BZQ','BZB','BZG','BZU', 'Beamed AGN', 'FSRQ', 'BZQ/Lense']

Fermi_Pulsars = ['PWN', 'PSR', 'pwn', 'psr', 'MSP']
Swift_Pulsars = ['Pulsar']

Fermi_Seyferts = ['nlsy1', 'sey', 'NLSY1']
Swift_Seyferts = ['Sy1','Sy1.9','Sy1.5','Sy2','Sy1.2','Sy1.8']

Fermi_HMXB = ['HMB', 'hmb', 'BIN']
Swift_HMXB = ['HMXB', 'XRB']

Fermi_LMXB = ['LMB', 'lmb', 'BIN', 'MSP']
Swift_LMXB = ['LMXB']

Fermi_Radio_Galaxy = ['RDG', 'rdg']
Swift_Radio_Galaxy = []

Fermi_Globular_Cluster = ['glc']
Swift_Globular_Cluster = ['GC']

Fermi_Galactic_Center = ['GC']
Swift_Galactic_Center = ['Galactic Center']

Fermi_Starburst_Galaxy = ['sbg']
Swift_Starburst_Galaxy = ['Starburst galaxy']

Fermi_SNR = ['snr']
Swift_SNR = ['SNR']

Fermi_Other_AGN = ['AGN', 'agn', 'css']
Swift_Other_AGN = ['Other AGN']

Fermi_Star = []
Swift_Star = ['CV']

Fermi_Types = [Fermi_Blazars, Fermi_Pulsars, Fermi_Seyferts, Fermi_HMXB, Fermi_LMXB, Fermi_Radio_Galaxy, Fermi_Globular_Cluster, Fermi_Galactic_Center, Fermi_Starburst_Galaxy, Fermi_SNR, Fermi_Other_AGN, Fermi_Star]
Swift_Types = [Swift_Blazars, Swift_Pulsars, Swift_Seyferts, Swift_HMXB, Swift_LMXB, Swift_Radio_Galaxy, Swift_Globular_Cluster, Swift_Galactic_Center, Swift_Starburst_Galaxy, Swift_SNR, Swift_Other_AGN, Swift_Star]

Fermi_Unknown = ['unk', '']
Swift_Unknown = ['Unknown', 'confused source', 'Unknown Star']

Flags = []

for i in range(len(All_Matches)):
    
    if All_Matches[i] in Both_Matched + Only_id_Matched:
        for x in range(len(Fermi_Types)):
            if fermi_source_types[i] in Fermi_Types[x] and swift_source_types[i] in Swift_Types[x]:
                Flags.append('M') 
            elif fermi_source_types[i] not in Fermi_Types[x] and swift_source_types[i] in Swift_Types[x]:
                if fermi_source_types[i] == 'agn' and swift_source_types[i] == 'Sy1':
                    Flags.append('M')
                elif fermi_source_types[i] not in Fermi_Unknown and swift_source_types[i] not in Swift_Unknown:
                    Flags.append('D')
                    
        if fermi_source_types[i] in Fermi_Unknown or swift_source_types[i] in Swift_Unknown:
            Flags.append('U')
                
    elif All_Matches[i] in Only_spatial_Matched:
        for x in range(len(Fermi_Types)):
            if fermi_source_types[i] in Fermi_Types[x] and swift_source_types[i] in Swift_Types[x]:
                Flags.append('M')
            elif fermi_source_types[i] not in Fermi_Types[x] and swift_source_types[i] in Swift_Types[x]:
                if fermi_source_types[i] in ['bll', 'bcu'] and swift_source_types[i] in ['Sy1', 'Sy1.9', 'Sy1.5']:
                    Flags.append('A')
                elif fermi_source_types[i] not in Fermi_Unknown and swift_source_types[i] not in Swift_Unknown:
                    Flags.append('F')
                    
        if fermi_source_types[i] in Fermi_Unknown or swift_source_types[i] in Swift_Unknown:
            Flags.append('U')
            

print(Counter(Flags))

Swift_Names, Fermi_Names = zip(*All_Matches)

All_Matched_Swift_RA = Spatial_Matches['Swift_RA'].tolist() + Id_Only_Matches['Swift_RA'].tolist()
All_Matched_Swift_Dec = Spatial_Matches['Swift_DEC'].tolist() + Id_Only_Matches['Swift_DEC'].tolist()
All_Matched_Fermi_RA = Spatial_Matches['Fermi_RA'].tolist() + Id_Only_Matches['Fermi_RA'].tolist()
All_Matched_Fermi_Dec = Spatial_Matches['Fermi_DEC'].tolist() + Id_Only_Matches['Fermi_DEC'].tolist()
All_Matched_Sep = Spatial_Matches['Sep'].tolist() + Id_Only_Matches['Sep'].tolist()

Match_Flag_Table = Table([Flags, Swift_Names, Fermi_Names, swift_source_types, fermi_source_types, All_Matched_Swift_RA, All_Matched_Swift_Dec, All_Matched_Fermi_RA, All_Matched_Fermi_Dec, All_Matched_Sep], names=['Flag', 'Swift_Name', 'Fermi_Name', 'Swift_Type', 'Fermi_Type', 'Swift_RA', 'Swift_DEC', 'Fermi_RA', 'Fermi_DEC', 'Sep'])
#Match_Flag_Table.pprint_all()
Match_Flag_Table.write(output_dir + '/Update_Match_Flagging.fits', overwrite=True)


"""

LAT_name_Minami = Minami_Flags['fermi_name_counterpart'].tolist()
BAT_name_Minami = Minami_Flags['bat_name_counterpart'].tolist()
LAT_class_Minami = Minami_Flags['fermi_category_type'].tolist()
BAT_class_Minami = Minami_Flags['bat_category_type'].tolist()

Minami_Classes = list(zip(BAT_class_Minami, LAT_class_Minami))
Tim_Classes = list(zip(swift_source_types, fermi_source_types))

Minami_Matches = list(zip(BAT_name_Minami, LAT_name_Minami))
Minami_Flags = Minami_Flags['flag_estimated'].tolist()

Minami_Match_Flags = list(zip(Minami_Flags, Minami_Matches, Minami_Classes))
Tim_Match_Flags = list(zip(Flags, All_Matches, Tim_Classes))

Minami_Match_Flags = list(zip(Minami_Flags, Minami_Matches, Minami_Classes))
Tim_Match_Flags = list(zip(Flags, All_Matches, Tim_Classes))

# Build name -> flag lookups
minami_dict = {match: flag for flag, match, classes in Minami_Match_Flags}
tim_dict = {match: flag for flag, match, classes in Tim_Match_Flags}

# Build name -> class lookups
minami_class_dict = {match: classes for flag, match, classes in Minami_Match_Flags}
tim_class_dict = {match: classes for flag, match, classes in Tim_Match_Flags}

# Check for matches present in both
common_matches = set(minami_dict.keys()) & set(tim_dict.keys())
only_in_minami = set(minami_dict.keys()) - set(tim_dict.keys())
only_in_tim = set(tim_dict.keys()) - set(minami_dict.keys())

print(f"Matches in both lists: {len(common_matches)}")
print(f"Only in Minami: {len(only_in_minami), only_in_minami}")
print(f"Only in Tim: {len(only_in_tim), only_in_tim}")

# Compare flags for matches present in both
mismatches = []
for match in common_matches:
    if minami_dict[match] != tim_dict[match]:
        mismatches.append((
            match,
            minami_dict[match], tim_dict[match],
            minami_class_dict[match], tim_class_dict[match]
        ))

print(f"\n{len(mismatches)} flag mismatches found:")
for match, minami_flag, tim_flag, minami_class, tim_class in mismatches:
    print(f"  {match, tim_class}: Minami={minami_flag}, Tim={tim_flag} ")
    
    
"""