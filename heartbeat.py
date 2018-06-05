from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from datetime import date
from scipy import stats, signal
import functions as funct
import os
from datetime import date
import matplotlib.dates as mdates


# Load mask data

mask_file = '/Users/jmh2g09/Documents/PhD/Data/Seasonal/seasonal_mask.nc'
nc = Dataset(mask_file, 'r')
mask = nc.variables['seasonal_mask'][:]
nc.close()

# Load SSH Data
DOT = np.full((59, 361, 60), fill_value = np.NaN)
dates = []
it = 0
for year in ['2011', '2012', '2013', '2014', '2015']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        dates.append(date(int(year), int(month), 1))
        file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc'
        if os.path.isfile(file):
            # Load the lat and lon data for each file
            nc = Dataset(file, 'r')
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            DOT[:, :, it] = nc.variables['dynamic_ocean_topography_seasonal_offset'][:] * mask
            nc.close()
            it += 1
            
DOT[DOT > -2] = np.NaN

heartbeat = np.nanmean(np.nanmean(DOT, 1), 0)

pl.figure()
pl.plot(dates, heartbeat)
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seasonal/heartbeat.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

# Load MCA stuff
# Import file
MCA_file = '/Users/jmh2g09/Documents/MATLAB/MCA_ssh_stresscurl.nc'
nc = Dataset(MCA_file, 'r')
ssh_PC1 = -nc.variables['ssh_PC1'][:]
ssh_PC2 = -nc.variables['ssh_PC2'][:]
wind_PC1 = -nc.variables['stress_PC1'][:]
wind_PC2 = -nc.variables['stress_PC2'][:]
EXPVAR = nc.variables['EXPVAR'][:]
nc.close()

A = 3.5
pl.figure()
fig, ax = pl.subplots()
pl.plot(dates, signal.detrend(heartbeat)*40, label='Mean SSHA', color='k', linestyle='-', linewidth=2)
pl.plot(dates, ssh_PC1+A, label='SSH PC(1)', color='r', linestyle='-', linewidth=1)
pl.plot(dates, wind_PC1+A, label='OSC PC(1)', color='r', linestyle='--', linewidth=1)
pl.plot(dates, ssh_PC2-A, label='SSH PC(2)', color='b', linestyle='-', linewidth=1)
pl.plot(dates, wind_PC2-A, label='OSC PC(2)', color='b', linestyle='--', linewidth=1)

pl.plot(dates, (EXPVAR[0]*ssh_PC1+EXPVAR[1]*ssh_PC2)/(EXPVAR[0]+EXPVAR[1]), color='purple', linestyle='-', linewidth=1)


# Shrink current axis by 20%
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
#ax.legend(loc='upper left', bbox_to_anchor=(1, 1))


pl.plot([date(2011, 8, 1), date(2011, 8, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)
pl.plot([date(2011, 11, 1), date(2011, 11, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)

pl.plot([date(2012, 9, 1), date(2012, 9, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)
pl.plot([date(2012, 11, 1), date(2012, 11, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)

pl.plot([date(2013, 9, 1), date(2013, 9, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)
pl.plot([date(2013, 12, 1), date(2013, 12, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)

pl.plot([date(2014, 9, 1), date(2014, 9, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)
pl.plot([date(2014, 11, 1), date(2014, 11, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)

pl.plot([date(2015, 7, 1), date(2015, 7, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)
pl.plot([date(2015, 10, 1), date(2015, 10, 1)], [-6.5, 6.5], color='k', linestyle='--', linewidth=0.5)

pl.xlim(date(2011, 1, 1), date(2016, 1, 1))
ax.xaxis.set_major_locator(mdates.YearLocator())
fig.autofmt_xdate()

pl.yticks([-4, 0, 4], ['PC(2)', 'Mean SSHA', 'PC(1)'])
pl.ylim([-6.5, 6.5])
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seasonal/heartbeat_plus_PC.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()






