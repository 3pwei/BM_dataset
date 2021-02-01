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
        
    for n in range(len(fn_temp)):
        
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
                    
        vol_uncertainty = vol_org*0.0169
        mv_err = np.nanmean(vol_uncertainty)
        err_list.append(mv_err)
        
    print(folder,"OK",len(data_list))
    return(data_list, err_list)

## import EXPT data
EXPT1 = uncertainty_mean_list("../ORC_Output/EXPT1/stomate/")
EXPT2 = uncertainty_mean_list("../ORC_Output/EXPT2/stomate/")

## plot figure ##
dataset = [EXPT1, EXPT2]
cases = ['EXPT1 (Transient simulation, Cyclic climate 1979-1994)','EXPT2 (Control simulation)']

fig = plt.figure(figsize=(12,8))

colors = ['#d62728',
          '#ff7f0e',
          '#2ca02c',
          '#1f77b4'] 

TFR=[1954,1972,1990,2012]
DFR=[94, 181,170, 228]

plt.bar(TFR[0],DFR[0], label=str('Forest Inventory'), width=-5, align='edge', alpha=0.2, color='#66B2FF')
plt.bar(TFR[1],DFR[1], width=-5, align='edge', alpha=0.2, color='#66B2FF')
plt.bar(TFR[2],DFR[2], width=-5, align='edge', alpha=0.2, color='#66B2FF')
plt.bar(TFR[3],DFR[3], width=-5, align='edge', alpha=0.2, color='#66B2FF')


T = list(range(1904,1980))
y = np.array(dataset[0][0])
plt.plot(T[:-1], y[:-1], marker='o', label=str(cases[0]), color='black')
plt.legend(bbox_to_anchor=(1.06, 1), loc='upper left', borderaxespad=0.)
y_err = np.array(dataset[0][1])
plt.fill_between(T, y - y_err, y + y_err, color='Grey', alpha=0.2)   

T = list(range(1979,2018))
y = np.array(dataset[1][0])
plt.plot(T, y, marker='o', label=str(cases[1]), color=colors[0])
plt.legend(bbox_to_anchor=(1.06, 1), loc='upper left', borderaxespad=0.)
y_err = np.array(dataset[1][1])
plt.fill_between(T, y - y_err, y + y_err, color='Grey', alpha=0.2)


plt.ylim((80,300))
plt.ylabel('Aboveground Wood Volume (m$^3$ ha$^{-1}$)', fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('Year', fontsize=16)
plt.xticks(fontsize=16)
plt.xlim((1910,2017))

plt.legend( fontsize=16, loc=2)

plt.text(TFR[0], DFR[0]+0.05, '%.0f' % TFR[0], ha='center', va= 'bottom',fontsize=16)
plt.text(TFR[1], DFR[1]+0.05, '%.0f' % TFR[1], ha='center', va= 'bottom',fontsize=16)
plt.text(TFR[2], DFR[2]+0.05, '%.0f' % TFR[2], ha='center', va= 'bottom',fontsize=16)
plt.text(TFR[3], DFR[3]+0.05, '%.0f' % TFR[3], ha='center', va= 'bottom',fontsize=16)

plt.grid(True)

plt.savefig("Fig_1a.tif", dpi = 100)
plt.show()
plt.close()


