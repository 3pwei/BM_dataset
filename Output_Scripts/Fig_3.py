# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import shapefile as shp


def uncertainty_mean_list(folder):
    
    fn = os.listdir("../Landcover_Map/")
    fn.sort()
    fn_temp = os.listdir(folder)
    fn_temp.sort()
    data_list = []
    err_list = []
    
    if len(fn) == len(fn_temp):
        pass
    else:
        print("number of files are difference.")
        return
        
    for n in range(len(fn)):
        
        temp =  xr.open_dataset(folder+fn_temp[n]).WOOD_VOLPIX[0].values*10000
        mv = np.nanmean(temp)
        data_list.append(mv)
        
        all_frac = xr.open_dataset("../Landcover_Map/"+fn[n]).maxvegetfrac.values[0, :, 32:118, 55:96]
        tree_frac = all_frac[1] + all_frac[2] + all_frac[3]
        vol_org = temp/tree_frac
        
        for i in range(vol_org.shape[0]):
            for j in range(vol_org.shape[1]):
                if np.isinf(vol_org[i,j]):
                    vol_org[i,j]=np.nan
                    
        vol_uncertainty = vol_org*0.0169  # Refer Chen et al.(2019)
        mv_err = np.nanmean(vol_uncertainty)
        err_list.append(mv_err)
        
    print(folder,"OK",len(data_list))
    return(data_list, err_list)

# import EXPT data
EXPT3A = uncertainty_mean_list("../ORC_Output/EXPT3A/stomate/")
EXPT3B = uncertainty_mean_list("../ORC_Output/EXPT3B/stomate/")
EXPT3C = uncertainty_mean_list("../ORC_Output/EXPT3C/stomate/")
EXPT2 = uncertainty_mean_list("../ORC_Output/EXPT2/stomate/")


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

## plot fiure ##
dataset = [EXPT2, EXPT3A, EXPT3B, EXPT3C]
cases = ['EXPT2 (Control simulation)', 'EXPT3A (1979 CO2)', 'EXPT3B (1979 LC)', 'EXPT3C (No windthrow)']

colors = ['#d62728',
          '#ff7f0e',
          '#2ca02c',
          '#1f77b4'] 

T = list(range(1979,2018))

plt.figure(figsize=(12,16))
plt.subplot(221)
#ax = fig.add_axes([0.1, 0.1, 1.8, 1.5])


for i in range(len(cases)):
    
    y = np.array(dataset[i][0])
    plt.plot(T, y, marker='o', label=str(cases[i]), color=colors[i])
    plt.legend(bbox_to_anchor=(1.06, 1), loc='upper left', borderaxespad=0.)
    
    y_err = np.array(dataset[i][1])
    plt.fill_between(T, y - y_err, y + y_err, color='Grey', alpha=0.2)


plt.ylim((190,260))
plt.ylabel('Aboveground Wood Volume (m$^3$ ha$^{-1}$)', fontsize=16)

plt.yticks(fontsize=12)
plt.xlabel('Year', fontsize=16)

plt.xticks(list(range(1980,2018,7)),fontsize=12)
plt.xlim((1979,2017))

plt.legend( fontsize=11, loc=2)

plt.grid(True)


# import spatail distribution
EXPT2H_f = xr.open_dataset('../ORC_Output/EXPT2/stomate/exp.t2_20170101_20171231_1M_stomate_history.nc').WOOD_VOLPIX[0].values*10000
EXPT3A_f = xr.open_dataset('../ORC_Output/EXPT3A/stomate/exp.t3.A_20170101_20171231_1M_stomate_history.nc').WOOD_VOLPIX[0].values*10000
EXPT3B_f = xr.open_dataset('../ORC_Output/EXPT3/EXPT3B/stomate/exp.t3.B_20170101_20171231_1M_stomate_history.nc').WOOD_VOLPIX[0].values*10000
EXPT3C_f = xr.open_dataset('../ORC_Output/EXPT3/EXPT3C/stomate/exp.t3.C_20170101_20171231_1M_stomate_history.nc').WOOD_VOLPIX[0].values*10000

#ColorBar
cm = LinearSegmentedColormap.from_list('name', ['#DBB700', 'white', '#00550E'])

# EXPT3A - EXPT2
SA = EXPT3A_f - EXPT2H_f 

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
        
plt.imshow(SA, origin='lower', vmin = -30, vmax = 30, cmap=cm,extent=extent)
cb1 = plt.colorbar(extend='both')
cb1.ax.tick_params(labelsize=12)
plt.axis('off')
plt.title('EXPT3A - Control Run',y=0,horizontalalignment='center', fontsize=16)

# EXPT3B - EXPT2
SB = EXPT3B_f - EXPT2H_f 

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
        
plt.imshow(SB, origin='lower', vmin = -300, vmax = 300, cmap=cm, extent=extent)
cb2 = plt.colorbar(extend='both')
cb2.ax.tick_params(labelsize=12)
plt.axis('off')
plt.title('EXPT3B - Control Run',y=0,horizontalalignment='center', fontsize=16)

# EXPT3C - EXPT2
SC = EXPT3C_f - EXPT2H_f 

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
        
cm = LinearSegmentedColormap.from_list('name', ['white', '#00550E'])
plt.imshow(SC, origin='lower', vmin = 0, vmax = 300, cmap=cm, extent=extent)
plt.clim(0, 300)
cb3 = plt.colorbar(extend='max')
cb3.ax.tick_params(labelsize=12)
plt.axis('off')
plt.title('EXPT3C - Control Run',y=0,horizontalalignment='center', fontsize=16)

plt.subplots_adjust(wspace =0, hspace =0.2)
plt.savefig(data+"Fig_3.tif", dpi = 100)
plt.show()
plt.close()
