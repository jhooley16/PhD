import os
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

year = '2015'
month = '03'

ssha = []
lon = []
lat = []
time = []
surface_type = []
retracking_mode = []

for mode in ['LRM', 'SAR', 'SARIN']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/' + year + month + '_' + mode)

    for file in os.listdir():

        f = open(file, 'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            # If data point is from open ocean (1) or from a lead (2)
            if columns[0] == '1' or columns[0] == '2':
                # If data point is listed as 'valid'
                if columns[1] == '1':
                    time.append(float(columns[4]))
                    lat.append(float(columns[5]))
                    lon.append(float(columns[6]))
                    ssha.append(float(columns[7]) - float(columns[8]))
                    surface_type.append(float(columns[0]))
                    if mode == 'LRM':
                        retracking_mode.append(1)
                    elif mode == 'SAR':
                        retracking_mode.append(2)
                    elif mode == 'SARIN':
                        retracking_mode.append(3)

ssha_sorted = [ssha for (time, ssha) in sorted(zip(time, ssha))]
lat_sorted = [lat for (time, lat) in sorted(zip(time, lat))]
lon_sorted = [lon for (time, lon) in sorted(zip(time, lon))]
surface_type_sorted = [surface_type for (time, surface_type) in sorted(zip(time, surface_type))]
retracking_mode_sorted = [retracking_mode for (time, retracking_mode) in sorted(zip(time, retracking_mode))]

pl.figure(figsize=(10, 10))
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
stereo_x, stereo_y = m(lon_sorted, lat_sorted)
m.scatter(stereo_x, stereo_y, c=retracking_mode_sorted, marker='.', edgecolors='none')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/'+ year + month + '_mode_coverage.png', format='png')
pl.close()

# Find the boundaries
LRM_2_SAR =[1] * 10 + [2] * 10
SAR_2_LRM = [2] * 10 + [1] * 10
LRM_2_SARIN = [1] * 10 + [3] * 10
SARIN_2_LRM = [3] * 10 + [1] * 10
SAR_2_SARIN = [2] * 10 + [3] * 10
SARIN_2_SAR = [3] * 10 + [2] * 10

iLRM2SAR = []
iLRM2SAR_lat = []
iLRM2SAR_lon = []
iSAR2LRM = []
iSAR2LRM_lat = []
iSAR2LRM_lon = []
iLRM2SARIN = []
iLRM2SARIN_lat = []
iLRM2SARIN_lon = []
iSARIN2LRM = []
iSARIN2LRM_lat = []
iSARIN2LRM_lon = []
iSAR2SARIN = []
iSAR2SARIN_lat = []
iSAR2SARIN_lon = []
iSARIN2SAR = []
iSARIN2SAR_lat = []
iSARIN2SAR_lon = []

for it in range(len(retracking_mode_sorted)):
    # LRM (1) to SAR (2)
    if retracking_mode_sorted[it:it + 20] == LRM_2_SAR:
        iLRM2SAR.append(it + 10)
        iLRM2SAR_lat.append(lat_sorted[it + 10])
        iLRM2SAR_lon.append(lon_sorted[it + 10])
    # SAR (2) to LRM (1)
    if retracking_mode_sorted[it:it + 20] == SAR_2_LRM:
        iSAR2LRM.append(it + 10)
        iSAR2LRM_lat.append(lat_sorted[it + 10])
        iSAR2LRM_lon.append(lon_sorted[it + 10])
    # LRM (1) to SARIN (3)
    if retracking_mode_sorted[it:it + 20] == LRM_2_SARIN:
        iLRM2SARIN.append(it + 10)
        iLRM2SARIN_lat.append(lat_sorted[it + 10])
        iLRM2SARIN_lon.append(lon_sorted[it + 10])
    # SARIN (3) to LRM (1)
    if retracking_mode_sorted[it:it + 20] == SARIN_2_LRM:
        iSARIN2LRM.append(it + 10)
        iSARIN2LRM_lat.append(lat_sorted[it + 10])
        iSARIN2LRM_lon.append(lon_sorted[it + 10])
    # SAR (2) to SARIN (3)
    if retracking_mode_sorted[it:it + 20] == SAR_2_SARIN:
        iSAR2SARIN.append(it + 10)
        iSAR2SARIN_lat.append(lat_sorted[it + 10])
        iSAR2SARIN_lon.append(lon_sorted[it + 10])
    # SARIN (3) to SAR (2)
    if retracking_mode_sorted[it:it + 20] == SARIN_2_SAR:
        iSARIN2SAR.append(it + 10)
        iSARIN2SAR_lat.append(lat_sorted[it + 10])
        iSARIN2SAR_lon.append(lon_sorted[it + 10])


pl.figure(figsize=(10, 10))
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])

# SAR-LRM boundaries
stereo_x, stereo_y = m(iLRM2SAR_lon, iLRM2SAR_lat)
m.scatter(stereo_x, stereo_y, marker='.', edgecolors='none', color='b')
stereo_x, stereo_y = m(iSAR2LRM_lon, iSAR2LRM_lat)
m.scatter(stereo_x, stereo_y, marker='.', edgecolors='none', color='b')

# SAR-SARIN boundaries
stereo_x, stereo_y = m(iSAR2SARIN_lon, iSAR2SARIN_lat)
m.scatter(stereo_x, stereo_y, marker='.', edgecolors='none', color='g')
stereo_x, stereo_y = m(iSARIN2SAR_lon, iSARIN2SAR_lat)
m.scatter(stereo_x, stereo_y, marker='.', edgecolors='none', color='g')

# LRM-SARIN boundaries
stereo_x, stereo_y = m(iLRM2SARIN_lon, iLRM2SARIN_lat)
m.scatter(stereo_x, stereo_y, marker='.', edgecolors='none', color='r')
stereo_x, stereo_y = m(iSARIN2LRM_lon, iSARIN2LRM_lat)
m.scatter(stereo_x, stereo_y, marker='.', edgecolors='none', color='r')

pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/'+ year + month + '_boundaries.png', format='png')
pl.close()


LRM2SAR_WEDD = []
LRM2SARIN_WEDD = []
SAR2SARIN_WEDD = []
LRM2SAR_IND = []
LRM2SARIN_IND = []
SAR2SARIN_IND = []
LRM2SAR_ROSS = []
LRM2SARIN_ROSS = []
SAR2SARIN_ROSS = []
LRM2SAR_AMBEL = []
LRM2SARIN_AMBEL = []
SAR2SARIN_AMBEL = []

for boundary in [iLRM2SAR, iSAR2LRM, iLRM2SARIN, iSARIN2LRM, iSAR2SARIN, iSARIN2SAR]:
    for step in boundary:
        if -60. <= np.mean(lon_sorted[step - 10:step + 10]) <= 0.:
            #print('Weddell')
            if boundary == iLRM2SAR:
                LRM2SAR_WEDD.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SAR_WEDD.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSAR2LRM:
                LRM2SAR_WEDD.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SAR_WEDD.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))

            if boundary == iLRM2SARIN:
                LRM2SARIN_WEDD.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SARIN_WEDD.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2LRM:
                LRM2SARIN_WEDD.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SARIN_WEDD.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
            
            if boundary==iSAR2SARIN:
                SAR2SARIN_WEDD.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_SAR2SARIN_WEDD.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2SAR:
                SAR2SARIN_WEDD.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_SAR2SARIN_WEDD.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))

        if 0 <= np.mean(lon_sorted[step - 10:step + 10]) <= 160.:
            #print('Indian')
            if boundary == iLRM2SAR:
                LRM2SAR_IND.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SAR_IND.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSAR2LRM:
                LRM2SAR_IND.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SAR_IND.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))

            if boundary == iLRM2SARIN:
                LRM2SARIN_IND.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SARIN_IND.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2LRM:
                LRM2SARIN_IND.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SARIN_IND.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
            
            if boundary==iSAR2SARIN:
                SAR2SARIN_IND.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_SAR2SARIN_IND.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2SAR:
                SAR2SARIN_IND.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_SAR2SARIN_IND.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
        
        if 160. <= np.mean(lon_sorted[step - 10:step + 10]) <= 180.:
            #print('Ross')
            if boundary == iLRM2SAR:
                LRM2SAR_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SAR_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSAR2LRM:
                LRM2SAR_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SAR_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))

            if boundary == iLRM2SARIN:
                LRM2SARIN_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SARIN_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2LRM:
                LRM2SARIN_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SARIN_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
            
            if boundary==iSAR2SARIN:
                SAR2SARIN_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_SAR2SARIN_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2SAR:
                SAR2SARIN_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_SAR2SARIN_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))

        if -180. <= np.mean(lon_sorted[step - 10:step + 10]) <= -130.:
            #print('Ross')
            if boundary == iLRM2SAR:
                LRM2SAR_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SAR_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSAR2LRM:
                LRM2SAR_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SAR_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))

            if boundary == iLRM2SARIN:
                LRM2SARIN_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SARIN_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2LRM:
                LRM2SARIN_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SARIN_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
            
            if boundary==iSAR2SARIN:
                SAR2SARIN_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_SAR2SARIN_ROSS.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2SAR:
                SAR2SARIN_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_SAR2SARIN_ROSS.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))

        if -130. <= np.mean(lon_sorted[step - 10:step + 10]) <= -60.:
            #print('Amundsen-Bellingshausen')
            if boundary == iLRM2SAR:
                LRM2SAR_AMBEL.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SAR_AMBEL.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSAR2LRM:
                LRM2SAR_AMBEL.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SAR_AMBEL.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))

            if boundary == iLRM2SARIN:
                LRM2SARIN_AMBEL.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_LRM2SARIN_AMBEL.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2LRM:
                LRM2SARIN_AMBEL.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_LRM2SARIN_AMBEL.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
            
            if boundary==iSAR2SARIN:
                SAR2SARIN_AMBEL.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
                #hist_SAR2SARIN_AMBEL.append(np.mean(ssha_sorted[step - 10:step]) - np.mean(ssha_sorted[step:step + 10]))
            
            if boundary==iSARIN2SAR:
                SAR2SARIN_AMBEL.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
                #hist_SAR2SARIN_AMBEL.append(np.mean(ssha_sorted[step:step + 10]) - np.mean(ssha_sorted[step - 10:step]))
