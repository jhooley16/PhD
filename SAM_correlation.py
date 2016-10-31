import numpy as np
import matplotlib.pyplot as pl
import os
from netCDF4 import Dataset
from scipy import stats
import functions as funct
from mpl_toolkits.basemap import Basemap

## Correlate the DOT data with the SAM index. 

# Load the SAM index file

sam_file = '/Users/jmh2g09/Documents/PhD/Data/Wind/SAM_2010-2016.txt'

sam_index = []
f = open(sam_file, 'r')
for line in f:
    line = line.strip()
    sam_index.append(float(line))
f.close()

sam_index = np.array(sam_index)

pl.figure()
pl.plot(sam_index)
pl.ylabel('SAM index')
pl.xlabel('Month (from NOV 2010)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Wind/Figures/SAM_index.png', transparent=True, dpi=300)
pl.close()

# Load the DOT data and calculate the correlation

## Open each file and append the DOT anomaly to a numpy array
# Initiate array of size of the data (59, 361) and the number of months (64)
dot_anom_ts = np.zeros((59, 361, len(sam_index)))
dot_2_anom_ts = np.zeros((59, 361, len(sam_index)))
# Initiate the index for the months
t = 0
months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

# Cycle through the years
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    # Cycle through the months
    for month in months:
        # If the file exists, open the file and extract the data
        file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year + '/Anomalies/' + year + month + '_DOT_anomaly.nc'
        if os.path.isfile(file):
            nc = Dataset(file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            # Assign the monthly gridded DOT anomaly to the array
            dot_anom_ts[:, :, t]  = nc.variables['dynamic_ocean_topography_anomaly'][:]
            dot_2_anom_ts[:, :, t]  = nc.variables['dynamic_ocean_topography_anomaly_no_offset'][:]
            nc.close()
            # Move the month index along by 1
            t += 1

# Calculate the correlation for each cell timeseries with the SAM index

## Cycle through each grid cell and calculate correlation
dot_anom_xcorr = np.full((59, 361), fill_value=np.NaN)
dot_anom_xcorr_pvalues = np.full((59, 361), fill_value=np.NaN)
dot_2_anom_xcorr = np.full((59, 361), fill_value=np.NaN)
dot_2_anom_xcorr_pvalues = np.full((59, 361), fill_value=np.NaN)

for ilat in range(len(lat)):
    for ilon in range(len(lon)):
        # If there are nans in the data
        if np.sum(np.isfinite(dot_anom_ts[ilat, ilon, :])) > len(dot_anom_ts[ilat, ilon, :]) // 2:
            xcorr = stats.spearmanr(dot_anom_ts[ilat, ilon, :], sam_index, nan_policy='omit')
            dot_anom_xcorr[ilat, ilon] = xcorr[0]
            dot_anom_xcorr_pvalues[ilat, ilon] = xcorr[1]
            xcorr_2 = stats.spearmanr(dot_2_anom_ts[ilat, ilon, :], sam_index, nan_policy='omit')
            dot_2_anom_xcorr[ilat, ilon] = xcorr_2[0]
            dot_2_anom_xcorr_pvalues[ilat, ilon] = xcorr_2[1]

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
grid_lats, grid_lons = np.meshgrid(lat, lon)
stereo_x, stereo_y = m(grid_lons, grid_lats)
        
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anom_xcorr)), cmap='RdBu_r')
m.colorbar()
pl.clim(1, -1)
m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anom_xcorr_pvalues)), [0.1], color='k')

#pl.title('Correlation between SAM and CryoSat-2 Altimetry')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Wind/Figures/SAM_correlation.png', transparent=True, dpi=300)
pl.close()

pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        
grid_lats, grid_lons = np.meshgrid(lat, lon)
stereo_x, stereo_y = m(grid_lons, grid_lats)
        
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_2_anom_xcorr)), cmap='RdBu_r')
m.colorbar()
pl.clim(1, -1)
m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_2_anom_xcorr_pvalues)), [0.1], color='k')

#pl.title('Correlation between SAM and CryoSat-2 Altimetry')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Wind/Figures/SAM_correlation_no_offset.png', transparent=True, dpi=300)
pl.close()