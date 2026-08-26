from astropy.table import Table
from matplotlib.colors import LogNorm
from astropy.coordinates import SkyCoord
from matplotlib.lines import Line2D

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u

output_dir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Outputs/'
catalog_dir = '/home/alab_student/Tim/Projects/Cross_Match_Update/Catalogs/'

cross_match_updated = Table.read(output_dir + 'cross_match_updated.fits')

fig = plt.figure(figsize=(15,9), facecolor='white')
ax = fig.add_subplot(111, projection="mollweide", facecolor='white')

ax.plot(np.zeros(100), np.linspace(-np.pi/2, np.pi/2, 100),
    color='black', linewidth=0.5)

ax.plot(np.linspace(-np.pi, np.pi, 500), np.zeros(500),
    color='black', linewidth=0.5)

Flags = ['M', 'D', 'F', 'A', 'U']
Labels = ['Firm Match (M)', 'Different Source Type (D)', 'False Match (F)', 'Ambiguous Association (A)', 'Unidentified Association (U)']
Marker = ['*', 'D', 'X', 'o', 'v']
Color = ['Blue', 'Green', 'Black', 'Red', 'Grey']
Sizes = [150,60,50,60,60]

legend_elements = [plt.Line2D([0], [0], marker=m, color='w', markerfacecolor=c, markeredgecolor=c, markersize=10, label=l) for m, c, l in zip(Marker, Color, Labels)]

for x in range(len(Flags)):
    mask = (cross_match_updated['flag'] == Flags[x])
    
    Swift_RA = cross_match_updated[mask]['bat_ra'].tolist()
    Swift_DEC = cross_match_updated[mask]['bat_dec'].tolist()
    
    RA = np.array(Swift_RA)
    DEC = np.array(Swift_DEC)

    coords_icrs = SkyCoord(ra=RA*u.deg, dec=DEC*u.deg, frame='icrs')
    coords_gal = coords_icrs.galactic

    l = coords_gal.l.wrap_at(180*u.deg).radian
    l=-l
    b = coords_gal.b.radian
    
    ax.scatter(l, b, marker=Marker[x], color=Color[x], s=Sizes[x], label=Labels[x], )

    for b_deg in [-90, -75, -60, -45, -30,-15, 0,15, 30,45, 60,75,90]:  
        b = np.deg2rad(b_deg)
        l_line = np.linspace(-np.pi, np.pi, 500)
        b_line = np.full_like(l_line, b)
        ax.plot(l_line, b_line, color='black', alpha=0.4, linewidth=1, linestyle=':')
        
    b_deg_list = [-75, -60, -45, -30,-15, 15, 30,45, 60,75]   
    b_rad_list = np.deg2rad(b_deg_list)

    for b_rad, b_deg in zip(b_rad_list, b_deg_list):
        ax.text(0.0, b_rad+0.07, f"{b_deg}°", color='black',
        ha='right', va='center', fontsize=12, zorder=20)
        
    for l_deg in [-150, -120, -90, -60, -30, 30, 60, 90, 120, 150]:
        l = -np.deg2rad(l_deg) 
        b_line = np.linspace(-np.pi/2, np.pi/2, 500)
        l_line = np.full_like(b_line, l)
        ax.plot(l_line, b_line, color='black', alpha=0.4, linewidth=1,linestyle=':')
        
    l_deg_list = [-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150]
    l_rad_list = -np.deg2rad(l_deg_list)  
    for l_rad, l_deg in zip(l_rad_list, l_deg_list):
            ax.text(l_rad - 0.07, 0.0, f"{l_deg}°", color='black',
            ha='center', va='bottom', fontsize=12, zorder=20)
        

    ax.set_xlabel(r"$l$ (°)")
    ax.set_ylabel(r"$b$ (°)")
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.1), ncol=5, fontsize=10, frameon=False, handletextpad=0.01, columnspacing=1.5)
    
plt.title('Update MeV Cross-Match All-sky Map', fontsize=15)

plt.tight_layout()

plt.savefig(output_dir + 'MeV_Updated_All_Sky_Map.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()