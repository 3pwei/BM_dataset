# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import shapefile as shp
from matplotlib import gridspec

def uncertainty_mean_list(folder):
    fn = os.listdir("/lfs/home/whuang/Landcover_Map/")
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
        
    for n in range(len(fn_temp)):
        
        temp =  xr.open_dataset(folder+fn_temp[n]).WOOD_VOLPIX[0].values*10000
        mv = np.nanmean(temp)
        data_list.append(mv)
        
        all_frac = xr.open_dataset("/lfs/home/whuang/Landcover_Map/"+fn[n]).maxvegetfrac.values[0, :, 32:118, 55:96]
        tree_frac = all_frac[1] + all_frac[2] + all_frac[3]
        vol_org = temp/tree_frac
        
        for i in range(vol_org.shape[0]):
            for j in range(vol_org.shape[1]):
                if np.isinf(vol_org[i,j]):
                    vol_org[i,j]=np.nan
                    
        vol_uncertainty = vol_org*0.0169
        mv_err = np.nanmean(vol_uncertainty)
        err_list.append(mv_err)
        
    print(folder,"OK",len(data_list))
    return(data_list, err_list)

## import shapefile ##
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

sf = shp.Reader(src+'SHP/Taiwan_part_84.shp')
sf2 = shp.Reader(src+'SHP/臺灣國家公園_wgs84.shp')

lower_x = 120.000244-0.05
lower_y = 21.520363
upper_x = 122.063416
upper_y = 25.489298
xedges = drange(lower_x, upper_x, 0.01) # 0.01 is the size of bin
yedges = drange(lower_y, upper_y, 0.01)
extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]


# import EXPT1 data
EXPT1 = uncertainty_mean_list("/lfs/home/whuang/RE_EXPT/EXPT1/stomate/")
# Remove outliers
def remove_outliers(input_array):
    input_array = input_array[input_array!=0]
    input_array = input_array[~np.isnan(input_array)]
    input_array = input_array[~np.isinf(input_array)]
    Q1,Q2,Q3 = np.percentile(input_array, [25, 50, 75])
    IQR = Q3-Q1
    input_array_max = Q3+1.5*IQR
    input_array_min = Q1-1.5*IQR
    return(float(input_array_max), float(input_array_min))

# import 2010 ABG data
temp =  xr.open_dataset('/lfs/home/whuang/RE_EXPT/EXPT1/stomate/re.exp.t1_19100101_19101231_1M_stomate_history.nc').WOOD_VOLPIX[0].values*10000
all_frac = xr.open_dataset("/lfs/home/whuang/Landcover_Map/PFT_map_tw_5km_1910.nc").maxvegetfrac.values[0, :, 32:118, 55:96]
tree_frac = all_frac[1] + all_frac[2] + all_frac[3]
vol_org = temp/tree_frac
for i in range(vol_org.shape[0]):
    for j in range(vol_org.shape[1]):
        if np.isinf(vol_org[i,j]) or vol_org[i,j]<0:
            vol_org[i,j]=np.nan

# import 1979 AGB data          
temp2 =  xr.open_dataset('/lfs/home/whuang/RE_EXPT/EXPT1/stomate/re.exp.t1_19790101_19791231_1M_stomate_history.nc').WOOD_VOLPIX[0].values*10000
all_frac2 = xr.open_dataset("/lfs/home/whuang/Landcover_Map/PFT_map_tw_5km_1979.nc").maxvegetfrac.values[0, :, 32:118, 55:96]
tree_frac2 = all_frac2[1] + all_frac2[2] + all_frac2[3]
vol_org2 = temp2/tree_frac2

## plot figure
dataset = [EXPT1]
cases = ['EXPT1 (Transient simulation)']

T = list(range(1904,1980))


fig = plt.figure(figsize=(24, 8)) 
gs = gridspec.GridSpec(1, 3, figure=fig) 
gs2 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0])
fig.add_subplot(gs2[:-1])

colors = ['#d62728',
          '#ff7f0e',
          '#2ca02c',
          '#1f77b4'] 
TFR=[1954,1972]
DFR=[94, 181]

plt.bar(TFR[0],DFR[0], label=str('Forest Inventory'), width=-5, align='edge', alpha=0.2, color='#66B2FF')
plt.bar(TFR[1],DFR[1], width=-5, align='edge', alpha=0.2, color='#66B2FF')

plt.text(TFR[0], DFR[0]+0.05, '%.0f' % TFR[0], ha='center', va= 'bottom',fontsize=16)
plt.text(TFR[1], DFR[1]+0.05, '%.0f' % TFR[1], ha='center', va= 'bottom',fontsize=16)

for i in range(len(cases)):
    
    y = np.array(dataset[i][0])
    
    plt.plot(T[6:], y[6:], marker='o', label=str(cases[i]), color='black')
    plt.legend(bbox_to_anchor=(1.06, 1), loc='upper left', borderaxespad=0.)
    
    y_err = np.array(dataset[i][1])
    plt.fill_between(T[6:], y[6:] - y_err[6:], y[6:] + y_err[6:], color='Grey', alpha=0.2)


plt.ylim((80,280))
plt.ylabel('Aboveground Wood Volume (m$^3$ ha$^{-1}$)', fontsize=16)
plt.yticks(fontsize=12)
plt.xlabel('Year', fontsize=16)
plt.xticks(fontsize=12)
plt.xlim((1910,1980))

plt.legend( fontsize=16, loc=2)
plt.grid(True)

gs4 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs2[2])
# Colorbar
levels = [ 0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0]
brewer_cmap = LinearSegmentedColormap.from_list('name', ['#FFFFFF', '#CCCC00','#00CC00', '#3333FF'], N=13)

fig.add_subplot(gs4[1])
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
plt.imshow(vol_org,origin='lower',vmax=650,vmin=0,cmap=brewer_cmap, extent=extent)
cb2=plt.colorbar(extend='max')
cb2.ax.tick_params(labelsize=12)
plt.axis('off')
plt.title('1910',y=-0.1,horizontalalignment='center', fontsize=16)

fig.add_subplot(gs4[0])
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

plt.imshow(vol_org2,origin='lower',vmax=650,vmin=0,cmap=brewer_cmap, extent=extent)
cb3=plt.colorbar(extend='max')
cb3.ax.tick_params(labelsize=12)
plt.axis('off')
plt.title('1979',y=-0.1,horizontalalignment='center', fontsize=16)
plt.tight_layout()
plt.savefig("Fig_S4.tif", dpi = 100)
plt.show()
plt.close()
