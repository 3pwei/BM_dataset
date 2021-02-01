# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt


def rainf_mean(path):
    fn = os.listdir(path)
    fn.sort()
    tw_map = np.zeros((86,41))
    for i in range(len(fn)):
        temp = xr.open_dataset(path+fn[i]).rainf[0]*365
        tw_map = tw_map + temp
    tw_map = tw_map/i
    return tw_map
def tair_mean(path):
    fn = os.listdir(path)
    fn.sort()
    tw_map = np.zeros((86,41))
    for i in range(len(fn)):
        temp = xr.open_dataset(path+fn[i]).tair[0]
        tw_map = tw_map + temp
    tw_map = tw_map/i
    return tw_map


## import EXPT data in 2009
rainf = rainf_mean('../ORC_Output/EXPT2/sechiba/')
tair = tair_mean('../ORC_Output/EXPT2/sechiba/')
WOOD_VOLPIX = xr.open_dataset('../ORC_Output/EXPT2/stomate/re.exp.t2_20090101_20091231_1M_stomate_history.nc').WOOD_VOLPIX[0]*10000
WOOD_VOL = xr.open_dataset('../ORC_Output/EXPT2/stomate/re.exp.t2_20090101_20091231_1M_stomate_history.nc').WOOD_VOL[0]*10000
PFT_2009 =xr.open_dataset('../ORC_Output/EXPT2/sechiba/re.exp.t2_20090101_20091231_1M_sechiba_history.nc').maxvegetfrac[0,:].values
WOOD_VOL = WOOD_VOL*PFT_2009

## import the 4th investigation of foreset
# Load 4th Investigation Data
BM_4th = pd.read_excel('../Investigate_Data/VOL_DBH_BN.xls')
temp = np.array((BM_4th.Ave_x.tolist(),BM_4th.Ave_y.tolist(),BM_4th.Ave_VOL_HA.tolist()))

TW_Map = np.zeros((86,41))
TW_Map[:]=np.nan
count = 0
for i in range(41):
    for j in range(86):
        
        jj=j+1
        ii=i+1
        if count<858:
            if temp.T[count][0]==ii and temp.T[count][1]==jj:
                TW_Map[j,i]=temp.T[count][2]
                count += 1
                
## processing the data
# Filter Mesk
Filter_Mesk1 = np.zeros((86,41))
Filter_Mesk1[:]=np.nan

for j in range(86):
    for i in range(41):
        if WOOD_VOLPIX[j,i].values>1:
            Filter_Mesk1[j,i]=1

Filter_Mesk2 = np.zeros((86,41))
Filter_Mesk2[:]=np.nan

for j in range(86):
    for i in range(41):
        if TW_Map[j,i]>0:
            Filter_Mesk2[j,i]=1    

PFT_2009 =xr.open_dataset('../ORC_Output/EXPT2/sechiba/re.exp.t2_20090101_20091231_1M_sechiba_history.nc').maxvegetfrac[0,:].values
# PFT Mesk
PFT_B = np.zeros((86,41))
PFT_N = np.zeros((86,41))
PFT_2 = np.zeros((86,41))
PFT_3 = np.zeros((86,41))
PFT_B[:] = np.nan
PFT_N[:] = np.nan
PFT_2[:] = np.nan
PFT_3[:] = np.nan


for j in range(86):
    for i in range(41):
        if (PFT_2009[1][j,i]>0.4 or PFT_2009[2][j,i]>0.4):
            PFT_B[j,i] = 1
        if (PFT_2009[3][j,i])>0.6:
            PFT_N[j,i] = 1
        if (PFT_2009[1][j,i])>0.4:
            PFT_2[j,i] = 1
        if (PFT_2009[2][j,i])>0.4:
            PFT_3[j,i] = 1   

# For Relationship between Precipitation and Boardleaf
tair_B_list = []
rainf_B_list = []
Investigation_B_list = []
ORC_B_list = []
ORC_VOL_B_list = []

#Filtering
tair_B = tair*Filter_Mesk1*Filter_Mesk2*PFT_B
rainf_B = rainf*Filter_Mesk1*Filter_Mesk2*PFT_B
TW_Map_B = TW_Map*Filter_Mesk1*Filter_Mesk2*PFT_B
WOOD_VOLPIX_B = WOOD_VOLPIX*Filter_Mesk1*Filter_Mesk2*PFT_B
WOOD_VOL_B = (WOOD_VOL[1]+WOOD_VOL[2])*Filter_Mesk1*Filter_Mesk2*PFT_B

for i in tair_B.values.tolist():
    for j in i:
        tair_B_list.append(j)
tair_B_list = np.array(tair_B_list)
tair_B_list = tair_B_list[~np.isnan(tair_B_list)]

for i in rainf_B.values.tolist():
    for j in i:
        rainf_B_list.append(j)
rainf_B_list = np.array(rainf_B_list)
rainf_B_list = rainf_B_list[~np.isnan(rainf_B_list)]

for i in TW_Map_B.tolist():
    for j in i:
        Investigation_B_list.append(j)
Investigation_B_list = np.array(Investigation_B_list)
Investigation_B_list = Investigation_B_list[~np.isnan(Investigation_B_list)]
        
for i in WOOD_VOLPIX_B.values.tolist():
    for j in i:
        ORC_B_list.append(j)
ORC_B_list = np.array(ORC_B_list)
ORC_B_list = ORC_B_list[~np.isnan(ORC_B_list)]

for i in WOOD_VOL_B.values.tolist():
    for j in i:
        ORC_VOL_B_list.append(j)
ORC_VOL_B_list = np.array(ORC_VOL_B_list)
ORC_VOL_B_list = ORC_VOL_B_list[~np.isnan(ORC_VOL_B_list)]

# For Relationship between Precipitation and Needleleaf
tair_N_list = []
rainf_N_list = []
Investigation_N_list = []
ORC_N_list = []
ORC_VOL_N_list = []

#Filtering
tair_N = tair*Filter_Mesk1*Filter_Mesk2*PFT_N
rainf_N = rainf*Filter_Mesk1*Filter_Mesk2*PFT_N
TW_Map_N = TW_Map*Filter_Mesk1*Filter_Mesk2*PFT_N
WOOD_VOLPIX_N = WOOD_VOLPIX*Filter_Mesk1*Filter_Mesk2*PFT_N
WOOD_VOL_N = WOOD_VOL[3]*Filter_Mesk1*Filter_Mesk2*PFT_N

for i in tair_N.values.tolist():
    for j in i:
        tair_N_list.append(j)
tair_N_list = np.array(tair_N_list)
tair_N_list = tair_N_list[~np.isnan(tair_N_list)]

for i in rainf_N.values.tolist():
    for j in i:
        rainf_N_list.append(j)
rainf_N_list = np.array(rainf_N_list)
rainf_N_list = rainf_N_list[~np.isnan(rainf_N_list)]

for i in TW_Map_N.tolist():
    for j in i:
        Investigation_N_list.append(j)
Investigation_N_list = np.array(Investigation_N_list)
Investigation_N_list = Investigation_N_list[~np.isnan(Investigation_N_list)]
        
for i in WOOD_VOLPIX_N.values.tolist():
    for j in i:
        ORC_N_list.append(j)
ORC_N_list = np.array(ORC_N_list)
ORC_N_list = ORC_N_list[~np.isnan(ORC_N_list)]

for i in WOOD_VOL_N.values.tolist():
    for j in i:
        ORC_VOL_N_list.append(j)
ORC_VOL_N_list = np.array(ORC_VOL_N_list)
ORC_VOL_N_list = ORC_VOL_N_list[~np.isnan(ORC_VOL_N_list)]

## plot figures ##
plt.figure(figsize=(11,8))

plt.subplot(221)
X  = rainf_B_list
Y1 = Investigation_B_list
Y2 = ORC_VOL_B_list
LHC_Rainf = rainf[51,17]
LHC_WOOD_VOL = 281.90
LHC_WOOD_VOL_ORC = WOOD_VOL[1,51,17]+WOOD_VOL[2,51,17]+WOOD_VOL[3,51,17]
HL_Rainf = rainf[45,28]
HL_WOOD_VOL = 100.38
HL_WOOD_VOL_ORC = WOOD_VOL[1,45,28]+WOOD_VOL[2,45,28]+WOOD_VOL[3,45,28]

plt.plot(X, Y1, "o", markersize=3, 
         markeredgewidth=1, markeredgecolor="#000000", markerfacecolor="#000000")
plt.plot(X, Y2, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="#A0A0A0", markerfacecolor="None")


plt.text(6500,930,r"$\bf{Broad leaf}$", ha='center')
plt.xlim(500,7500)
plt.ylim(0,1000)
plt.xlabel('Annual Precipitation (mm)')
plt.ylabel('Aboveground Wood Volume ($m^3$ ha$^{-1}$)')

plt.twinx()
plt.plot(LHC_Rainf, LHC_WOOD_VOL, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="r")
plt.plot(HL_Rainf, HL_WOOD_VOL, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="r")
plt.plot(LHC_Rainf, LHC_WOOD_VOL_ORC, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="None")
plt.plot(HL_Rainf, HL_WOOD_VOL_ORC, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="None")
plt.plot([LHC_Rainf,LHC_Rainf], [LHC_WOOD_VOL+8,LHC_WOOD_VOL_ORC-8], 'r--',linewidth=0.5)
plt.plot([HL_Rainf,HL_Rainf], [HL_WOOD_VOL+8,HL_WOOD_VOL_ORC-8], 'r--',linewidth=0.5)

plt.ylim(0,1000)
plt.yticks([LHC_WOOD_VOL,HL_WOOD_VOL], ['Site 1','Site 2'])




plt.subplot(222)
X  = rainf_N_list
Y1 = Investigation_N_list
Y2 = ORC_VOL_N_list
CLM_Rainf = rainf[65,27]
CLM_WOOD_VOL = 275.29
CLM_WOOD_VOL_ORC = WOOD_VOL[1,65,27]+WOOD_VOL[2,65,27]+WOOD_VOL[3,65,27]

plt.plot(X, Y1, "o", markersize=3, 
         markeredgewidth=1, markeredgecolor="#000000", markerfacecolor="#000000")
plt.plot(X, Y2, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="#A0A0A0", markerfacecolor="None")

plt.text(6500,930,r"$\bf{Needle leaf}$", ha='center')
plt.xlim(500,7500)
plt.ylim(0,1000)
plt.xlabel('Annual Precipitation (mm)')
plt.ylabel('Aboveground Wood Volume ($m^3$ ha$^{-1}$)')

plt.twinx()
plt.plot(CLM_Rainf, CLM_WOOD_VOL, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="r")
plt.plot(CLM_Rainf, CLM_WOOD_VOL_ORC, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="None")
plt.plot([CLM_Rainf,CLM_Rainf], [CLM_WOOD_VOL+8,CLM_WOOD_VOL_ORC-8], 'r--',linewidth=0.5)


plt.ylim(0,1000)
plt.yticks([CLM_WOOD_VOL], ['Site 3'])



plt.subplot(223)
X  = tair_B_list
Y1 = Investigation_B_list
Y2 = ORC_VOL_B_list
LHC_Tair = tair[51,17]
LHC_WOOD_VOL = 281.90
LHC_WOOD_VOL_ORC = WOOD_VOL[1,51,17]+WOOD_VOL[2,51,17]+WOOD_VOL[3,51,17]
HL_Tair = tair[45,28]
HL_WOOD_VOL = 100.38
HL_WOOD_VOL_ORC = WOOD_VOL[1,45,28]+WOOD_VOL[2,45,28]+WOOD_VOL[3,45,28]

plt.plot(X, Y1, "o", markersize=3, 
         markeredgewidth=1, markeredgecolor="#000000", markerfacecolor="#000000")
plt.plot(X, Y2, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="#A0A0A0", markerfacecolor="None")

plt.xlim(286,306)
plt.ylim(0,1000)
plt.xlabel('Annual Temperature (K)')
plt.ylabel('Aboveground Wood Volume ($m^3$ ha$^{-1}$)')

plt.twinx()
plt.plot(LHC_Tair, LHC_WOOD_VOL, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="r")
plt.plot(HL_Tair, HL_WOOD_VOL, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="r")
plt.plot(LHC_Tair, LHC_WOOD_VOL_ORC, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="None")
plt.plot(HL_Tair, HL_WOOD_VOL_ORC, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="None")
plt.plot([LHC_Tair,LHC_Tair], [LHC_WOOD_VOL+8,LHC_WOOD_VOL_ORC-8], 'r--',linewidth=0.5)
plt.plot([HL_Tair,HL_Tair], [HL_WOOD_VOL+8,HL_WOOD_VOL_ORC-8], 'r--',linewidth=0.5)

plt.ylim(0,1000)
plt.yticks([LHC_WOOD_VOL,HL_WOOD_VOL], ['Site 1','Site 2'])


plt.subplot(224)
X  = tair_N_list
Y1 = Investigation_N_list
Y2 = ORC_VOL_N_list
CLM_Tair = tair[65,27]
CLM_WOOD_VOL = 275.29
CLM_WOOD_VOL_ORC = WOOD_VOL[1,65,27]+WOOD_VOL[2,65,27]+WOOD_VOL[3,65,27]


plt.plot(X, Y1, "o", markersize=3, 
         markeredgewidth=1, markeredgecolor="#000000", markerfacecolor="#000000")
plt.plot(X, Y2, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="#A0A0A0", markerfacecolor="None")

#plt.legend(['Observation', 'ORCHIDEE'], title=r"$\bf{Needle leaf}$", framealpha=1)
plt.xlim(286,306)
plt.ylim(0,1000)
plt.xlabel('Annual Temperature (K)')
plt.ylabel('Aboveground Wood Volume ($m^3$ ha$^{-1}$)')

plt.twinx()
plt.plot(CLM_Tair, CLM_WOOD_VOL, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="r")
plt.plot(CLM_Tair, CLM_WOOD_VOL_ORC, "o", markersize=4, 
         markeredgewidth=0.8, markeredgecolor="r", markerfacecolor="None")
plt.plot([CLM_Tair,CLM_Tair], [CLM_WOOD_VOL+8,CLM_WOOD_VOL_ORC-8], 'r--',linewidth=0.5)

plt.ylim(0,1000)
plt.yticks([CLM_WOOD_VOL], ['Site 3'])



plt.tight_layout(pad=3.0) # subplot distance
plt.savefig("Fig_2b.tiff", dpi = 100)



plt.show()     
