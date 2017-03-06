from netCDF4 import Dataset
import numpy as np
from scipy import interpolate, signal, stats
from datetime import date
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import os
import functions as funct

## Load GRACE data
GRACE = np.full((64, 59, 361), fill_value=np.NaN)
it = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        GRACE_file = '/Users/jmh2g09/Documents/PhD/Data/GRACE/' + year + month + '_GRACE.nc'
        if os.path.isfile(GRACE_file):
            nc = Dataset(GRACE_file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            GRACE[it, :, :] = nc.variables['GRACE'][:]
            nc.close()
            it += 1

for it in range(64):
    GRACE[it, :, :] = GRACE[it, :, :]

# Load the SAM index file
sam_file = '/Users/jmh2g09/Documents/PhD/Data/Wind/SAM_2010-2016.txt'
sam_index = []
f = open(sam_file, 'r')
for line in f:
    line = line.strip()
    sam_index = np.append(sam_index, float(line))
f.close()

## Correlation between GRACE and SAM

GRACE_xcorr = np.full((59, 361), fill_value=np.NaN)
GRACE_xcorr_pvalues = np.full((59, 361), fill_value=np.NaN)

for ilat in range(len(lat)):
    for ilon in range(len(lon)):
        if np.sum(np.isfinite(GRACE[:, ilat, ilon])) > len(GRACE[:, ilat, ilon]) // 2:
            ts_1 = GRACE[:, ilat, ilon]
            xcorr = stats.pearsonr(ts_1, sam_index)
            GRACE_xcorr[ilat, ilon] = xcorr[0]
            GRACE_xcorr_pvalues[ilat, ilon] = xcorr[1]

#GRACE_xcorr[57, 300] = 1000

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
        
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(GRACE_xcorr)), cmap='RdBu_r')
m.colorbar()
pl.clim(1, -1)
m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(GRACE_xcorr_pvalues)), [0.2], color='k')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/GRACE/Validation/GRACE_SAM_corr.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()