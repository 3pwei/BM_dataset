# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import shapefile as shp
from matplotlib.colors import LinearSegmentedColormap

# import land-cover maps
PFT_1980 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_1980.nc').maxvegetfrac[0,1:4,32:118, 55:96].values
PFT_1980 = np.nansum(PFT_1980,axis=0)
PFT_1995 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_1995.nc').maxvegetfrac[0,1:4,32:118, 55:96].values
PFT_1995 = np.nansum(PFT_1995,axis=0)
PFT_2000 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_2000.nc').maxvegetfrac[0,1:4,32:118, 55:96].values
PFT_2000 = np.nansum(PFT_2000,axis=0)
PFT_2005 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_2005.nc').maxvegetfrac[0,1:4,32:118, 55:96].values
PFT_2005 = np.nansum(PFT_2005,axis=0)
PFT_2010 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_2010.nc').maxvegetfrac[0,1:4,32:118, 55:96].values
PFT_2010 = np.nansum(PFT_2010,axis=0)
PFT_2015 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_2015.nc').maxvegetfrac[0,1:4,32:118, 55:96].values
PFT_2015 = np.nansum(PFT_2015,axis=0)


# plot the figure 7
# setting the figure size
plt.figure(figsize=(12,16))

cm = LinearSegmentedColormap.from_list('name', ['#DBB700', 'white', '#00550E'])

plt.subplot(221)
# laod Taiwan shapefile
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
        
BN_9580 = PFT_1995-PFT_1980
plt.imshow(BN_9580,origin='lower',cmap=cm,vmax=1,vmin=-1, extent=extent)
plt.colorbar(shrink=0.8)
plt.axis('off')
plt.title('1980 to 1995', fontsize=16)

plt.subplot(222)
# laod Taiwan shapefile
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
        

BN_2095 = PFT_2000-PFT_1995
plt.imshow(BN_2095,origin='lower',cmap=cm,vmax=1,vmin=-1, extent=extent)
plt.colorbar(shrink=0.8)
plt.axis('off')
plt.title('1995 to 2000', fontsize=16)

plt.subplot(223)
# laod Taiwan shapefile
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
        
BN_0500 = PFT_2005-PFT_2000
plt.imshow(BN_0500,origin='lower',cmap=cm,vmax=0.4,vmin=-0.4, extent=extent)
plt.colorbar(shrink=0.8)
plt.axis('off')
plt.title('2000 to 2005', fontsize=16)

plt.subplot(224)
# laod Taiwan shapefile
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
        
BN_1005 = PFT_2010-PFT_2005
plt.imshow(BN_1005,origin='lower',cmap=cm,vmax=0.4,vmin=-0.4, extent=extent)
plt.colorbar(shrink=0.8)
plt.axis('off')
plt.title('2005 to 2010', fontsize=16)

plt.savefig("Fig_7.tif", dpi = 100)
