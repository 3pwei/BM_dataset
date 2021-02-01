# This script was written by Wei Huang
# Contact email: popop9987@gmail.com

import numpy as np
import xarray as xr
import os
import pandas as pd
import matplotlib.pyplot as plt

# dtm
dataset_2d = np.array(pd.read_excel('../DTM/WRF_DTM_5km_T97.xlsx')[['y','x','Ave_GRIDCO']])
taiwan_2d = np.array(pd.read_excel('../DTM/WRF_taiwan_raster_5km_T97.xlsx')[['y','x']])

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


# Import EXPT2 data
EXPT2_2017 = xr.open_dataset('../EXPT2/stomate/re.exp.t2_20170101_20171231_1M_stomate_history.nc')
EXPT2_2017_ABG = EXPT2_2017.WOOD_VOLPIX[0].values*10000

# Define a tool for the pie plot
def drawPieMarker(ratios, xs, ys, sizes, colors):
    assert sum(ratios) <= 1, 'sum of ratios needs to be < 1'

    markers = []
    previous = 0
    # calculate the points of the pie pieces
    for color, ratio in zip(colors, ratios):
        this = 2 * np.pi * ratio + previous
        x  = [0] + np.cos(np.linspace(previous, this, 10)).tolist() + [0]
        y  = [0] + np.sin(np.linspace(previous, this, 10)).tolist() + [0]
        xy = np.column_stack([x, y])
        previous = this
        markers.append({'marker':xy, 's':np.abs(xy).max()**2*np.array(sizes), 'facecolor':color})

    # scatter each of the pie pieces to create pies
    for marker in markers:
        ax.scatter(xs, ys, **marker)


def draw_pie(dist, 
             xpos, 
             ypos, 
             size, 
             ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))

    # for incremental pie slices
    cumsum = np.cumsum(dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()

    for r1, r2 in zip(pie[:-1], pie[1:]):
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.cos(angles).tolist()
        y = [0] + np.sin(angles).tolist()

        xy = np.column_stack([x, y])
        ax.scatter([xpos], [ypos], marker=xy, s=size, c=["r", "b", "g"])

    return ax

# plot figures
x = plt.subplots(figsize=(8,8))

#composition
COM = EXPT2_2017.WOOD_VOL[0].values
COM_SUM = COM[1]+COM[2]+COM[3]
BF = (COM[1]+COM[2])/COM_SUM
NF = COM[3]/COM_SUM

dtm_list = []
for i in range(dtm_array.shape[0]):
    for j in range(dtm_array.shape[1]):
        dtm_list.append(dtm_array[i,j])
        
ABG_list = []
for i in range(dtm_array.shape[0]):
    for j in range(dtm_array.shape[1]):
        ABG_list.append(EXPT2_2017_ABG[i,j])
BF_list = []
for i in range(dtm_array.shape[0]):
    for j in range(dtm_array.shape[1]):
        BF_list.append(BF[i,j])    
NF_list = []
for i in range(dtm_array.shape[0]):
    for j in range(dtm_array.shape[1]):
        NF_list.append(NF[i,j])    

for h in range(len(ABG_list)):
    if BF_list[h] == BF_list[h]:
        drawPieMarker(ratios=[BF_list[h],NF_list[h]], xs = dtm_list[h], ys = ABG_list[h], sizes = 100, colors = ['g','y'])


plt.xlabel('Elevation (m, a.s.l.)',fontsize=16)
plt.xticks(fontsize=12)
plt.ylabel('AGB (m$^3$ ha$^{-1}$)',fontsize=16)
plt.yticks(fontsize=12)
plt.grid(True, linestyle='--')
plt.title('Control Run',fontsize=16)
plt.savefig("Fig_8.tif", dpi = 100)
plt.show()

