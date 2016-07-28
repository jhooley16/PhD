import numpy as np
import os
import matplotlib.pyplot as pl
from netCDF4 import Dataset
import functions as funct
from mpl_toolkits.basemap import Basemap
from scipy import stats

mode = input('What mode? (full, summer, or winter): ')

if mode == 'full':
    time = np.arange(1, 65) / 12
    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
elif mode == 'summer':
    time = np.array((1, 2, 3, 4, 13, 14, 15, 16, 25, 26, 27, 28, 37, 38, 39, 40, 49, 50, 51, 52, 61, 62)) / 12
    months = ['01', '02', '03', '04']
elif mode == 'winter':
    time = np.array((1, 2, 7, 8, 9, 10, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36, 37, 38, 43, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 58, 59, 60, 61, 62)) / 12
    months = ['05', '06', '07', '08', '09', '10', '11', '12']

## Open each file and append the DOT anomaly to a numpy array
# Initiate array of size of the data (59, 361) and the number of months (64)
ssh_anom_ts = np.zeros((59, 361, len(time)))
# Initiate the index for the months
t = 0

# Cycle through the years
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    # Change the directory
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/SSH/' + year + '/Anomalies')
    # Cycle through the months
    for month in months:
        # If the file exists, open the file and extract the data
        file = year + month + '_SSH_anomaly.nc'
        if os.path.isfile(file):
            nc = Dataset(file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            # Assign the monthly gridded DOT anomaly to the array
            ssh_anom_ts[:, :, t]  = nc.variables['sea_surface_height_anomaly'][:]
            nc.close()
            # Move the month index along by 1
            t += 1

## Cycle through each grid cell and calculate the trend of the data
ssh_anom_regress = np.full((59, 361), fill_value=np.NaN)

for ilat in range(len(lat)):
    for ilon in range(len(lon)):
        # If there are nans in the data
        if np.any(np.isfinite(ssh_anom_ts[ilat, ilon, :])):
            regress = stats.linregress(time, funct.inpaint_nans(ssh_anom_ts[ilat, ilon, :]))
            ssh_anom_regress[ilat, ilon] = (regress[0] * 1000)

## Calculate the mean circumpolar trend

# Calculate the grid size for each cell
grid_lon, grid_lat = np.meshgrid(lon, lat)
# Calculate the surface area of each cell
S = funct.surface_area(grid_lat, grid_lon, 0.5, 1.0)

total_area = np.nansum(np.nansum(~np.isnan(ssh_anom_regress) * S))

total_trend = np.nansum(ssh_anom_regress * S) / total_area

print('The total trend is:', total_trend, 'mm / year')

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
#m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
grid_lats, grid_lons = np.meshgrid(lat, lon)
stereo_x, stereo_y = m(grid_lons, grid_lats)
        
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ssh_anom_regress)), cmap='RdBu_r')
m.colorbar()
pl.clim(20, -20)
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seasonal/' + mode + '_SSH_trend_map.png', transparent=True, dpi=300)
pl.close()