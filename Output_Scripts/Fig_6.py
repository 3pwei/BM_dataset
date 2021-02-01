# This script was written by Wei Huang
# Contact email: popop9987@gmail.com


import numpy as np
import xarray as xr
import os
import shapefile as shp
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

## import wind data
Wind_daily = xr.open_dataset('../ORC_output/EXPT2_daily/stomate/exp.t2.daily_20170101_20171231_1M_stomate_history.nc').DAILY_MAX_WIND

def windstorm_acc(array3D):
    windstorm_acc = np.zeros((86,41))
    windstorm_acc[:] = np.nan
    for j in range(array3D.shape[1]):
        for i in range(array3D.shape[2]):
            acc = array3D[:,j,i][array3D[:,j,i]>20].shape[0]
            windstorm_acc[j,i] = acc
    return windstorm_acc

def Total_winds_acc(Folder):
    fn = os.listdir(Folder)
    fn.sort()
    total_acc = []
    for i in range(len(fn)):
        Wind_daily = xr.open_dataset('../ORC_output/EXPT2_daily/stomate/'+fn[i]).DAILY_MAX_WIND
        total_acc.append(windstorm_acc(Wind_daily))
        print(fn[i])
    total_acc = np.array(total_acc)
    total_accf = np.nanmean(total_acc,axis=0)/365
    return total_accf

## processing the wind data
test = Total_winds_acc('../ORC_output/EXPT2_daily/stomate/')

## import shapefile
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

sf = shp.Reader('../SHP/Taiwan_part_84.shp')
sf2 = shp.Reader('../SHP/臺灣國家公園_wgs84.shp')

lower_x = 120.000244-0.05
lower_y = 21.520363
upper_x = 122.063416
upper_y = 25.489298

# histogram
xedges = drange(lower_x, upper_x, 0.01) # 0.01 is the size of bin
yedges = drange(lower_y, upper_y, 0.01)
extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]


## plot figure ##
plt.figure(figsize=(6,9))
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
cm = LinearSegmentedColormap.from_list('name', ['white', 'orange'])
plt.imshow(test, origin='lower', vmin = 0, vmax = 0.01, cmap=cm, extent=extent)
cbar = plt.colorbar( shrink=0.8,extend='max', format=OOMFormatter(-2, mathText=True))
cbar.ax.tick_params(labelsize=12)
cbar.update_ticks()
font = {'family' : 'normal',
        'color'  : 'k',
        'weight' : 'normal',
        'size'   : 12,
        }
cbar.set_label('Frequency (1x10$^{-2}$ day$^{-1}$)', rotation=90,fontdict=font, labelpad=10,fontsize=14)

plt.axis('off')
plt.title('MaxWindSpeed > 20 m/s',y=0,horizontalalignment='center', fontsize=16)
plt.savefig(data+"Fig_6.tiff", dpi = 100)
plt.show()
plt.close()
