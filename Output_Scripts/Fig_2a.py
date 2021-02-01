# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import numpy as np
import xarray as xr
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.cm as mpl_cm
import shapefile as shp
from matplotlib.colors import LinearSegmentedColormap


## Taiwan shapefile Tools
def drange(start, stop, step):
    result = list()
    r = start
    i = 0.0
    while r < stop:
        result.append(r)
        i += 1.0
        r = round(start + i * step, 3)
    result.append(r)
    return result

# Read shapefile
sf = shp.Reader('../SHP/Taiwan_part_84.shp')
sf2 = shp.Reader('../SHP/臺灣國家公園_wgs84.shp')

# Set boundaries
lower_x = 120.000244-0.05
lower_y = 21.520363
upper_x = 122.063416
upper_y = 25.489298

# Set variable for extent
xedges = drange(lower_x, upper_x, 0.01) # 0.01 is the size of bin
yedges = drange(lower_y, upper_y, 0.01)
extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]


EXPT2_2017 = xr.open_dataset('../EXPT2/sechiba/re.exp.t2_20170101_20171231_1M_sechiba_history.nc')
nav_lon = EXPT2_2017.nav_lon.values
nav_lat = EXPT2_2017.nav_lat.values

## plot figure ##
# setting figure size
fig = plt.figure(figsize=(8,8))

# laoding Taiwan shapefile
for shape in sf.shapeRecords():
    # end index of each components of map
    l = shape.shape.parts

    len_l = len(l)  # how many parts of countries i.e. land and islands
    x = [i[0] for i in shape.shape.points[:]] # list of latitude
    y = [i[1] for i in shape.shape.points[:]] # list of longitude
    l.append(len(x)) # ensure the closure of the last component
    for k in range(len_l):
        # draw each component of map.
        # l[k] to l[k + 1] is the range of points that make this component
        plt.plot(x[l[k]:l[k + 1]],y[l[k]:l[k + 1]], 'k-',linewidth=1)

for shape in sf2.shapeRecords():
    # end index of each components of map
    l = shape.shape.parts

    len_l = len(l)  # how many parts of countries i.e. land and islands
    x = [i[0] for i in shape.shape.points[:]] # list of latitude
    y = [i[1] for i in shape.shape.points[:]] # list of longitude
    l.append(len(x)) # ensure the closure of the last component
    for k in range(len_l):
        # draw each component of map.
        # l[k] to l[k + 1] is the range of points that make this component
        plt.plot(x[l[k]:l[k + 1]],y[l[k]:l[k + 1]], 'grey',linewidth=1)

# setting the colorbar levels and colors
levels = [ 0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0]
brewer_cmap = LinearSegmentedColormap.from_list('name', ['#FFFFFF', '#CCCC00','#00CC00', '#3333FF'], N=13)

# plotting the 2D figure
WVP = plt.imshow(park_mark, origin='lower', vmax=650, vmin=0, cmap=brewer_cmap, extent=extent)
 
# plotting dashed lines
plt.scatter(nav_lon[45,28],nav_lat[45,28], facecolors='none', color = ['k'])
plt.text((nav_lon[45,31]+nav_lon[44,31])/2,(nav_lat[45,31]+nav_lat[44,31])/2, "Site 1", color='k',fontsize=12)
plt.plot([nav_lon[45,28]+0.05,nav_lon[45,31]-0.02], [nav_lat[45,28],nav_lat[45,31]], 'k--',linewidth=1)

plt.scatter(nav_lon[51,17],nav_lat[51,17], facecolors='none', color = ['k'])
plt.text((nav_lon[51,4]+nav_lon[50,4])/2,(nav_lat[51,5]+nav_lat[50,5])/2, "Site 2", color='k',horizontalalignment='right',fontsize=12)
plt.plot([nav_lon[51,4]+0.02,nav_lon[51,17]-0.05], [nav_lat[51,4],nav_lat[51,17]], 'k--',linewidth=1)

plt.scatter(nav_lon[65,27],nav_lat[65,27], facecolors='none', color = ['k'])
plt.text((nav_lon[65,38]+nav_lon[64,38])/2,(nav_lat[65,38]+nav_lat[64,38])/2, "Site 3", color='k',fontsize=12)
plt.plot([nav_lon[65,27]+0.05,nav_lon[65,38]-0.02], [nav_lat[65,27],nav_lat[65,38]], 'k--',linewidth=1)

# adding city's name
plt.text(nav_lon[75,26],nav_lat[75,26], "Taipei", color='r',fontsize=12)
plt.text(nav_lon[58,11],nav_lat[58,11], "Taichung", color='r',fontsize=12)
plt.text(nav_lon[33,2],nav_lat[33,2], "Tainan", color='r',fontsize=12)
plt.axis('off')

# setting colorbar location
cbaxes = fig.add_axes([0.05, 0.18, 0.8, 0.4], frameon=False) 
cbar = plt.colorbar(WVP, shrink=0.9,extend='max')
font = {'family' : 'normal',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 12,
        }
cbar.set_label('Aboveground Wood Volume (m$^3$ ha$^{-1}$)', rotation=90,fontdict=font, labelpad=-55,fontsize=11)

plt.axis('off')
plt.savefig("Fig_2a.tif", dpi = 100)
plt.show()
plt.close()
