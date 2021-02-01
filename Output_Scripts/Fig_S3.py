# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import numpy as np
import xarray as xr
import pandas as pd
import shapefile as shp
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

# import DTM
dataset_2d = np.array(pd.read_excel('../DTM/WRF_DTM_5km_T97.xlsx')[['y','x','Ave_GRIDCO']])
taiwan_2d = np.array(pd.read_excel('../DTM/WRF_taiwan_raster_5km_T97.xlsx')[['y','x']])

# DTM spatial distribution
dtm_array = np.zeros((86,41))
dtm_array[:] = np.nan

for i in range(len(dataset_2d)):
    max_values = np.nanmax(dataset_2d.T[:][2])
    dtm_array[int(dataset_2d[i][0])-1, int(dataset_2d[i][1])-1] = dataset_2d[i][2]

tw_array = np.zeros((86,41))
tw_array[:] = np.nan

for j in range(len(taiwan_2d)):
    tw_array[int(taiwan_2d[j][0])-1, int(taiwan_2d[j][1])-1] = 1

dtm_array = dtm_array*tw_array

# remove outliers
def remove_outliers(input_array):
    input_array_A = input_array
    input_array = input_array[input_array!=0]
    input_array = input_array[~np.isnan(input_array)]
    input_array = input_array[~np.isinf(input_array)]
    Q1,Q2,Q3 = np.percentile(input_array, [25, 50, 75])
    IQR = Q3-Q1
    input_array_max = Q3+1.5*IQR
    input_array_min = Q1-1.5*IQR
    
    # rebuild array
    New_A = np.zeros((input_array_A.shape))
    New_A[:] = np.nan
    for j in range(input_array_A.shape[0]):
        for i in range(input_array_A.shape[1]):
            if input_array_A[j,i] > input_array_max:
                New_A[j,i] = 0
            elif input_array_A[j,i] < input_array_min:
                New_A[j,i] = 0
            else:
                New_A[j,i] = input_array_A[j,i]
    
    return(float(input_array_max), float(input_array_min), New_A)

## import data of temperature and precipitation
fn = os.listdir('../ORC_Output/EXPT2/sechiba/')
fn.sort()

Tair = []
Precip = []
for i in range(len(fn)):
    Tair.append(remove_outliers(xr.open_dataset('../ORC_Output/EXPT2/sechiba/'+fn[i]).tair[0].values)[2])
    Precip.append(remove_outliers(xr.open_dataset('../ORC_Output/EXPT2/sechiba/'+fn[i]).rainf[0].values*365)[2])
Tair = np.array(Tair)
Precip = np.array(Precip)
Tair_m = np.nanmean(Tair,axis=0)
Precip_m = np.nanmean(Precip,axis=0)

Tair_1979 = xr.open_dataset('../ORC_Output/EXPT2/sechiba/exp.t2_19790101_19791231_1M_sechiba_history.nc').tair[0].values
Tair_2017 = xr.open_dataset('../ORC_Output/EXPT2/sechiba/exp.t2_20170101_20171231_1M_sechiba_history.nc').tair[0].values
Precip_1979 = xr.open_dataset('../ORC_Output/EXPT2/sechiba/exp.t2_19790101_19791231_1M_sechiba_history.nc').rainf[0].values*365
Precip_2017 = xr.open_dataset('../ORC_Output/EXPT2/sechiba/exp.t2_20170101_20171231_1M_sechiba_history.nc').rainf[0].values*365
Tair_7917 = remove_outliers(Tair_1979)[2] - remove_outliers(Tair_2017)[2]
Precip_7917 = remove_outliers(Precip_1979)[2] - remove_outliers(Precip_2017)[2]

## import shapefiles
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

## plot figure ##
plt.figure(figsize=(16,8))

plt.subplot(131)
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
cm_DEM = LinearSegmentedColormap.from_list('name', ['white', 'black'])
plt.imshow(dtm_array,origin='lower', cmap=cm_DEM, vmax=3000, vmin=0, extent=extent)
plt.colorbar(extend='max', shrink=0.8)
plt.axis('off')
plt.title('Digital Elevation Model (m)', fontsize=16)


plt.subplot(132)
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
        
#ColorBar
cm_P = LinearSegmentedColormap.from_list('name', ['white','#66B2FF','#66FF66','yellow','orange', '#FF3333'])

plt.imshow(Tair_m,origin='lower', cmap=cm_P,vmax=305,vmin=275, extent=extent)
plt.colorbar(extend='max', shrink=0.8)
plt.axis('off')
plt.title('Temprature (K)', fontsize=16)


plt.subplot(133)
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
cm_P = LinearSegmentedColormap.from_list('name', ['white','#66B2FF','#66FF66','yellow','orange', '#FF3333'])
plt.imshow(Precip_m,origin='lower', cmap=cm_P,vmax=3000,vmin=500, extent=extent)
plt.colorbar(extend='max', shrink=0.8)
plt.axis('off')
plt.title('Aunnal precipitation (mm)', fontsize=16)
plt.savefig("Fig_S3.tif", dpi = 100)
plt.show()
