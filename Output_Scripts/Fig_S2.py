# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import shapefile as shp
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.basemap import Basemap


# import the land-cover  maps
PFT_1979 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_1979.nc')
PFT_2009 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_2009.nc')
PFT_2017 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_2017.nc')

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

sf = shp.Reader('/data1/home/whuang/i_python/BM/src/SHP/Taiwan_part_84.shp')
sf2 = shp.Reader('/data1/home/whuang/i_python/BM/src/SHP/臺灣國家公園_wgs84.shp')

lower_x = 120.000244-0.05
lower_y = 21.520363
upper_x = 122.063416
upper_y = 25.489298

# histogram
xedges = drange(lower_x, upper_x, 0.01) # 0.01 is the size of bin
yedges = drange(lower_y, upper_y, 0.01)
extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]

## plot the figure S2##
# setting the figure size
plt.figure(figsize=(8,16))
# setting the colorbar
cm = LinearSegmentedColormap.from_list('name', [ 'white', 'green'])

## plotting the sub-figure 1
plt.subplot(321)
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
        
B1979 = PFT_1979.maxvegetfrac[0,1,32:118, 55:96].values + PFT_1979.maxvegetfrac[0,2,32:118, 55:96].values
plt.imshow(B1979,origin='lower',cmap=cm,vmax=1,vmin=0, extent=extent)
plt.colorbar()
plt.axis('off')
plt.title('1979 BroadLeavedTropical')

## plotting the sub-figure 2
plt.subplot(322)
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
        
N1979 = PFT_1979.maxvegetfrac[0,3,32:118, 55:96].values
plt.imshow(N1979,origin='lower',cmap=cm,vmax=1,vmin=0, extent=extent)
plt.colorbar()
plt.axis('off')
plt.title('1979 NeedleleafTemperate')

## plotting the sub-figure 3
plt.subplot(323)
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
        
B2009 = PFT_2009.maxvegetfrac[0,1,32:118, 55:96].values + PFT_2009.maxvegetfrac[0,2,32:118, 55:96].values
plt.imshow(B2009,origin='lower',cmap=cm,vmax=1,vmin=0, extent=extent)
plt.colorbar()
plt.axis('off')
plt.title('2009 BroadLeavedTropical')

## plotting the sub-figure 4
plt.subplot(324)
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
        
N2009 = PFT_2009.maxvegetfrac[0,3,32:118, 55:96].values
plt.imshow(N2009,origin='lower',cmap=cm,vmax=1,vmin=0, extent=extent)
plt.colorbar()
plt.axis('off')
plt.title('2009 NeedleleafTemperate')

## plotting the sub-figure 5
plt.subplot(325)
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
        
B2017 = PFT_2017.maxvegetfrac[0,1,32:118, 55:96].values + PFT_2017.maxvegetfrac[0,2,32:118, 55:96].values
plt.imshow(B2017,origin='lower',cmap=cm,vmax=1,vmin=0, extent=extent)
plt.colorbar()
plt.axis('off')
plt.title('2017 BroadLeavedTropical')

## plotting the sub-figure 6
plt.subplot(326)
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
        
N2017 = PFT_2017.maxvegetfrac[0,3,32:118, 55:96].values
plt.imshow(N2017,origin='lower',cmap=cm,vmax=1,vmin=0, extent=extent)
plt.colorbar()
plt.axis('off')
plt.title('2017 NeedleleafTemperate')
plt.savefig("Fig_S2.tiff", dpi = 100)
