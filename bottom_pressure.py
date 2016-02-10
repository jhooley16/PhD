import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl

# Bottom pressure recorder located in the Weddell Sea
# -60.8249 N
# -54.7221 E (305.2779 W)
# depth = 1920 m 
file = '/Users/jmh2g09/Documents/PhD/Data/BPR/DPS_DEEP_1113_DQ105443_drp.txt'

day = []
bpa_notide_nodrift = []

f = open(file, 'r')

for line in f:
    line = line.strip()
    columns = line.split()
    if np.size(columns) > 4:
        # only choose data for the year 2012
        if columns[2] == '2012':
            # only choose 'valid' data
            if columns[1] == '0':
                day.append(float(columns[3]))
                bpa_notide_nodrift.append(float(columns[5]))

# bottom pressure is measured in millibars, need to convert to Pascals
# 1 millibar = 100 Pa
bpa_pascals = np.array(bpa_notide_nodrift) * 100

# Calculate the approximate change in sea level necessary to produce the pressure change
dh = - bpa_pascals / (1028 * 9.81)

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/2012/MDT/Anomalies')

mdta_at_recorder = []

for file in os.listdir():
    if file[-14:] == 'mdt_anomaly.nc':
        nc = Dataset(file, 'r')
        mdta = nc.variables['Mean Dynamic Height Anomaly'][:]
        lat = nc.variables['Latitude'][:]
        lon = nc.variables['Longitude'][:]
        nc.close()
        
        ilat = np.where(lat == -61)
        ilon = np.where(lon == -55)
        
        mdta_at_recorder.append(float(mdta[ilon, ilat]))

mdt_day = np.arange(15, 360, 30)

pl.figure()
pl.plot(day, dh)
pl.plot(mdt_day, np.array(mdta_at_recorder))
pl.plot()
pl.show()