# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt

# Remnve outliers
def remove_outliers(array):
    test=[]
    for i in array:
        for j in i:
            test.append(j)

    test = np.array(test)
    test = test[test!=0]
    test = test[~np.isnan(test)]
    test = test[~np.isinf(test)]
    Q1,Q2,Q3 = np.percentile(test, [25, 50, 75])
    IQR = Q3-Q1
    test_max = Q3+1.5*IQR
    test_min = Q1-1.5*IQR
    test1 = test[test>test_min]
    test2 = test1[test1<test_max]
    return(float(test_max), float(test_min),test2)

# uncertainty mean function
def WUE_list(folder):
    
    fn_temp_sechiba = os.listdir(folder+'/sechiba/')
    fn_temp_sechiba.sort()
    fn_temp_stomate = os.listdir(folder+'/stomate/')
    fn_temp_stomate.sort()
    fn_temp_PFT = os.listdir('../Landcover_Map/')
    fn_temp_PFT.sort()
    
    data_list = []
    
    if len(fn_temp_sechiba) == len(fn_temp_stomate)==39:
        pass
    else:
        print("number of files are difference.")
        return
        
    for n in range(39):
        PFT = xr.open_dataset('../Landcover_Map/'+fn_temp_PFT[n])
        sechiba =  xr.open_dataset(folder+'/sechiba/'+fn_temp_sechiba[n])
        stomate =  xr.open_dataset(folder+'/stomate/'+fn_temp_stomate[n])
        tran = (sechiba.transpir[0,1].values/PFT.maxvegetfrac[0,1,32:118, 55:96].values+
                sechiba.transpir[0,2].values/PFT.maxvegetfrac[0,2,32:118, 55:96].values+
                sechiba.transpir[0,3].values/PFT.maxvegetfrac[0,3,32:118, 55:96].values)/3
        GPP = (stomate.GPP[0,1].values+
               stomate.GPP[0,2].values+
               stomate.GPP[0,3].values)/3
        WUE = GPP/tran
        
        data_max,data_min,data_filter = remove_outliers(WUE)
        
        temp = data_filter
        mv = np.nanmean(temp)
        data_list.append(mv)
        
        
    print(folder,"OK",len(data_list))
    return(data_list)


# rainfall
def rain_mean_list(folder):
    fn_temp_sechiba = os.listdir(folder+'/sechiba/')
    fn_temp_sechiba.sort()
    fn_temp_stomate = os.listdir(folder+'/stomate/')
    fn_temp_stomate.sort()
    data_list = []

    for n in range(39):
        
        sechiba =  xr.open_dataset(folder+'/sechiba/'+fn_temp_sechiba[n])
        stomate =  xr.open_dataset(folder+'/stomate/'+fn_temp_stomate[n])
        rainf = sechiba.rainf[0]*365
        GPP = (stomate.GPP[0,1]+stomate.GPP[0,2]+stomate.GPP[0,3])/3
        
        # mesk for forsest land
        Mesk = np.zeros((86,41))
        Mesk[:]=np.nan
        for j in range(GPP.shape[0]):
            for i in range(GPP.shape[1]):
                if GPP[j,i]>0:
                    Mesk[j,i] = 1
  
        # Filter by Mesk
        rainf_m = rainf*Mesk
        
        temp = rainf_m
        mv = np.nanmean(temp)
        data_list.append(mv)
        
        
    print(folder,"OK",len(data_list))
    return(data_list)

## import EXPT data
EXPT3A = WUE_list("../ORC_Output/EXPT3A/")
EXPT3B = WUE_list("../ORC_Output/EXPT3B/")
EXPT3C = WUE_list("../ORC_Output/EXPT3C/")
EXPT2 = WUE_list("../ORC_Output/EXPT2/")
EXPT2_Rain = rain_mean_list("../ORC_Output/EXPT2/")

## plot figure ##
from sklearn.metrics import r2_score
import itertools
#sort data by precip
Data = np.array([EXPT2_Rain,EXPT3B])
Data2 = sorted(Data.T, key = lambda s: s[0])
Data2 = np.array(Data2)

x0=Data2.T[0]
y0=Data2.T[1]

for i in list(itertools.combinations(range(1,40),2)):
    index = list(i)
    x = np.delete(x0[2:],index)
    y = np.delete(y0[2:],index)
    z = np.polyfit(x, y, 1)
    y_hat = np.poly1d(z)(x)
    plt.plot(x, y_hat, "#C0C0C0", lw=0.5)
    #plt.gca().text(0.07, 0.93, "N = 2",transform=plt.gca().transAxes,
    #fontsize=12, verticalalignment='top')

z = np.polyfit(x0[2:], y0[2:], 1)
y_hat = np.poly1d(z)(x0[2:])
plt.plot(x0[2:], y_hat, "r--", lw=1)

plt.plot(x0[2:],y0[2:],"+", ms=8, mec="k", markerfacecolor="None")
plt.plot(x0[0:2],y0[0:2],"+", ms=8, mec="#808080", markerfacecolor="None")

plt.xlabel('Annual Precipitation (mm)')
plt.ylabel('Water Use Efficiency (g m$^{-2}$mm$^{-1}$)')    
plt.grid(True)

plt.savefig("Fig_5.tiff", dpi = 100)
plt.show()
