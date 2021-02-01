# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import shapefile as shp

## import EXPT data
sechiba_2 = xr.open_dataset('../ORC_Output/EXPT2/sechiba/re.exp.t2_20170101_20171231_1M_sechiba_history.nc')
stomate_2 = xr.open_dataset('../ORC_Output/EXPT2/stomate/re.exp.t2_20170101_20171231_1M_stomate_history.nc')

sechiba_A = xr.open_dataset('../ORC_Output/EXPT3A/sechiba/re.exp.t3.A_20170101_20171231_1M_sechiba_history.nc')
stomate_A = xr.open_dataset('../ORC_Output/EXPT3A/stomate/re.exp.t3.A_20170101_20171231_1M_stomate_history.nc')

sechiba_B = xr.open_dataset('../ORC_Output/EXPT3B/sechiba/re.exp.t3.B_20170101_20171231_1M_sechiba_history.nc')
stomate_B = xr.open_dataset('../ORC_Output/EXPT3B/stomate/re.exp.t3.B_20170101_20171231_1M_stomate_history.nc')

sechiba_C = xr.open_dataset('../ORC_Output/EXPT3C/sechiba/re.exp.t3.C_20170101_20171231_1M_sechiba_history.nc')
stomate_C = xr.open_dataset('../ORC_Output/EXPT3C/stomate/re.exp.t3.C_20170101_20171231_1M_stomate_history.nc')

PFT_2017 = xr.open_dataset('../Landcover_Map/PFT_map_tw_5km_2017.nc')

tran_2 = (sechiba_2.transpir[0,1].values/PFT_2017.maxvegetfrac[0,1,32:118, 55:96].values+
          sechiba_2.transpir[0,2].values/PFT_2017.maxvegetfrac[0,2,32:118, 55:96].values+
          sechiba_2.transpir[0,3].values/PFT_2017.maxvegetfrac[0,3,32:118, 55:96].values)/3
tran_A = (sechiba_A.transpir[0,1].values/PFT_2017.maxvegetfrac[0,1,32:118, 55:96].values+
          sechiba_A.transpir[0,2].values/PFT_2017.maxvegetfrac[0,2,32:118, 55:96].values+
          sechiba_A.transpir[0,3].values/PFT_2017.maxvegetfrac[0,3,32:118, 55:96].values)/3
tran_B = (sechiba_B.transpir[0,1].values/PFT_2017.maxvegetfrac[0,1,32:118, 55:96].values+
          sechiba_B.transpir[0,2].values/PFT_2017.maxvegetfrac[0,2,32:118, 55:96].values+
          sechiba_B.transpir[0,3].values/PFT_2017.maxvegetfrac[0,3,32:118, 55:96].values)/3
tran_C = (sechiba_C.transpir[0,1].values/PFT_2017.maxvegetfrac[0,1,32:118, 55:96].values+
          sechiba_C.transpir[0,2].values/PFT_2017.maxvegetfrac[0,2,32:118, 55:96].values+
          sechiba_C.transpir[0,3].values/PFT_2017.maxvegetfrac[0,3,32:118, 55:96].values)/3

GPP_2 = (stomate_2.GPP[0,1].values/PFT_2017.maxvegetfrac[0,1,32:118, 55:96].values+
         stomate_2.GPP[0,2].values/PFT_2017.maxvegetfrac[0,2,32:118, 55:96].values+
         stomate_2.GPP[0,3].values/PFT_2017.maxvegetfrac[0,3,32:118, 55:96].values)/3
WUE_2 = GPP_2/tran_2/10

GPP_A = (stomate_A.GPP[0,1].values/PFT_2017.maxvegetfrac[0,1,32:118, 55:96].values+
         stomate_A.GPP[0,2].values/PFT_2017.maxvegetfrac[0,2,32:118, 55:96].values+
         stomate_A.GPP[0,3].values/PFT_2017.maxvegetfrac[0,3,32:118, 55:96].values)/3
WUE_A = GPP_A/tran_A/10

GPP_B = (stomate_B.GPP[0,1].values/PFT_2017.maxvegetfrac[0,1,32:118, 55:96].values+
         stomate_B.GPP[0,2].values/PFT_2017.maxvegetfrac[0,2,32:118, 55:96].values+
         stomate_B.GPP[0,3].values/PFT_2017.maxvegetfrac[0,3,32:118, 55:96].values)/3
WUE_B = GPP_B/tran_B/10

GPP_C = (stomate_C.GPP[0,1].values/PFT_2017.maxvegetfrac[0,1,32:118, 55:96].values+
         stomate_C.GPP[0,2].values/PFT_2017.maxvegetfrac[0,2,32:118, 55:96].values+
         stomate_C.GPP[0,3].values/PFT_2017.maxvegetfrac[0,3,32:118, 55:96].values)/3
WUE_C = GPP_C/tran_C/10


# Remove outliers
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

# Proceesing the ORC outout data to the average values
EXPT3A = WUE_list("../ORC_Output/EXPT3A/")
EXPT3B = WUE_list("../ORC_Output/EXPT3B/")
EXPT3C = WUE_list("../ORC_Output/EXPT3C/")
EXPT2 = WUE_list("../ORC_Output/EXPT2/")
EXPT2_Rain = rain_mean_list("../ORC_Output/EXPT2/")

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



## plot figures
from matplotlib.colors import LinearSegmentedColormap

dataset = [EXPT2, EXPT3A, EXPT3B, EXPT3C]
cases = ['EXPT2 (Control simulation)', 'EXPT3A (1979 CO2)', 'EXPT3B (1979 LC)', 'EXPT3C (No windthrow)']

colors = ['#d62728',
          '#ff7f0e',
          '#2ca02c',
          '#1f77b4'] 

T = list(range(1979,2018))

plt.figure(figsize=(12,16))
ax1 = plt.subplot(221)

for i in range(len(cases)):
    
    y = np.array(dataset[i])
    plt.plot(T, y, marker='*', markersize=5, label=str(cases[i]), color=colors[i], zorder=2600)
    #ax.legend(bbox_to_anchor=(1.06, 1), loc=2, borderaxespad=0.)


plt.ylabel('Water Use Efficiency (g m$^{-2}$ mm$^{-1}$)', fontsize=16)
plt.yticks(fontsize=12)
plt.xlabel('Year', fontsize=16)
plt.ylim((2.0,5.))
plt.xticks(list(range(1980,2018,7)),fontsize=12)
plt.xlim((1979,2017))
#plt.legend(fontsize=10, loc=2)
plt.grid(True)

# twin object for two different y-axis on the sample plot
ax2 = ax1.twinx()
y2= np.array(EXPT2_Rain)
plt.bar(T, y2 , label='Precipitation', color='Gray', alpha=0.5 )
plt.ylabel('Annual Precipitation (mm)', fontsize=16, labelpad=-20)
plt.ylim((0,24000))
plt.yticks(list(range(0,6000,2000)),fontsize=12)
#plt.legend(fontsize=12, loc=3)

h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1+h2, l1+l2, loc=2,fontsize=14)

#ColorBar
#cm = LinearSegmentedColormap.from_list('name', ['blue', 'white', 'red'])
cm = LinearSegmentedColormap.from_list('name', ['#DBB700', 'white', '#00550E'])

# EXPT3A - EXPT2
SA = WUE_A - WUE_2

# Remove nan & inf & outliers
data_max, data_min, data_filter = remove_outliers(SA)
SA[SA>data_max]=0
SA[SA<data_min]=0

plt.subplot(222)
#cm = plt.cm.get_cmap('RdBu')
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
        
plt.imshow(SA, origin='lower', cmap=cm, vmax=3, vmin=-3,extent=extent)
#plt.imshow(SA, origin='lower', cmap=cm, vmax=3, vmin=-3)
cb1 = plt.colorbar(extend='both')
cb1.ax.tick_params(labelsize=12)
plt.axis('off')
plt.title('EXPT3A - Control Run',y=0,horizontalalignment='center', fontsize=16)


# EXPT3B - EXPT2
SB = WUE_B - WUE_2

# Remove nan & inf & outliers
data_max, data_min, data_filter = remove_outliers(SB)
SB[SB>data_max]=0
SB[SB<data_min]=0

plt.subplot(223)
#cm = plt.cm.get_cmap('RdBu')

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
        
plt.imshow(SB, origin='lower', cmap=cm, vmax=3, vmin=-3,extent=extent)
cb2 = plt.colorbar(extend='both')
cb2.ax.tick_params(labelsize=12)
plt.axis('off')
plt.title('EXPT3B - Control Run',y=0,horizontalalignment='center', fontsize=16)

# EXPT3C - EXPT2
SC = WUE_C - WUE_2

# Remove nan & inf & outliers
data_max, data_min, data_filter = remove_outliers(SC)
SC[SC>data_max]=0
SC[SC<data_min]=0

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
        
#cm = LinearSegmentedColormap.from_list('name', ['blue', 'white', 'red'])
#cm = LinearSegmentedColormap.from_list('name', ['white', '#00550E'])
plt.imshow(SC, origin='lower', cmap=cm,  vmax=3, vmin=-3, extent=extent)
cb3 = plt.colorbar(extend='both')
cb3.ax.tick_params(labelsize=12)
plt.axis('off')
plt.title('EXPT3C - Control Run',y=0,horizontalalignment='center', fontsize=16)

plt.subplots_adjust(wspace =0.2, hspace =0.2)

plt.savefig(data+"Fig_4.tif", dpi = 100)
plt.show()
plt.close()
