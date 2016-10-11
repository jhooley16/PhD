import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from scipy import stats
from mpl_toolkits.basemap import Basemap

## Create a bottom pressure/tide gauge - CS2 correlation map

# Load the bottom pressure data
for station_name  in  ['AntBasePrat_bpr', 'Argentine_Islands_tide_gauge',
    'Drake_passage_north_deep', 'Drake_passage_north', 'Drake_passage_south_deep',
    'Drake_passage_south', 'Syowa_station_bpr']:

    pressure_file = '/Users/jmh2g09/Documents/PhD/Data/BPR/Processed/' + station_name + '_processed.csv'

    year = []
    month = []
    bpr = []

    f = open(pressure_file, 'r')
    for line in f:
        line = line.strip()
        columns = line.split(',')
        year.append(int(columns[0]))
        month.append(int(columns[1]))
        bpr.append(float(columns[2]))
    f.close()

    locations_file = '/Users/jmh2g09/Documents/PhD/Data/BPR/Processed/locations.csv'

    # Open the locations (lat, lon) for the stations, to mark on correlation map
    f = open(locations_file, 'r')
    for line in f:
        line = line.strip()
        columns = line.split(',')
        if columns[0] == station_name:
            print(columns[0])
            station_lat = float(columns[1])
            station_lon = float(columns[2])
    f.close()

    # Load the altimetry data for the timeseries defined by the bottom
    # pressure/tide gauge record

    dot_anom_ts = np.zeros((59, 361, len(year)))

    for t in range(len(year)):

        # Make the month string consistent with the file convention
        if int(month[t]) < 10:
            month_str = '0' + str(month[t])
        elif int(month[t]) >= 10:
            month_str = str(month[t])

        altimetry_file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' \
            + str(year[t]) + '/Anomalies/' + str(year[t]) + month_str + '_DOT_anomaly.nc'

        nc = Dataset(altimetry_file, 'r')
        lat = nc.variables['latitude'][:]
        lon = nc.variables['longitude'][:]
        dot_anom_ts[:, :, t] = nc.variables['dynamic_ocean_topography_anomaly'][:]
        nc.close()

    # Calculate the correlation between the time series and the altimetry data
    dot_anom_xcorr = np.full((59, 361), fill_value=np.NaN)
    dot_anom_xcorr_pvalues = np.full((59, 361), fill_value=np.NaN)

    for ilat in range(len(lat)):
        for ilon in range(len(lon)):
            # If there are nans in the altimetry data
            if np.any(np.isfinite(dot_anom_ts[ilat, ilon, :])):
                xcorr = stats.spearmanr(funct.inpaint_nans(dot_anom_ts[ilat, ilon, :]), bpr)
                dot_anom_xcorr[ilat, ilon] = xcorr[0]
                dot_anom_xcorr_pvalues[ilat, ilon] = xcorr[1]

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

    location_x, location_y = m(station_lon, station_lat)
    m.scatter(location_x, location_y, marker='*', s=50, color='k', zorder=100)

    pl.title('CryoSat-2 Altimetry and ' + station_name.replace('_', ' ') + ' Correlation')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/BPR/Figures/' + station_name + '_correlation.png',
        transparent=True, dpi=300)
    pl.close()