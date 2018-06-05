import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from datetime import date
import matplotlib.dates as mdates
from scipy import stats, signal
from mpl_toolkits.basemap import Basemap

## Create a bottom pressure/tide gauge - CS2 correlation map

bathy_file = '/Users/jmh2g09/Documents/PhD/Data/Bathymetry/GEBCO_2014_2D.nc'
nc = Dataset(bathy_file, 'r')
bathy_lon = nc.variables['lon'][:]
bathy_lat = nc.variables['lat'][:]
bathy_data = nc.variables['elevation'][:]
nc.close()

nc = Dataset('INPUT_bathy.nc', 'w')
nc.createDimension('lat', np.size(bathy_lat))
nc.createDimension('lon', np.size(bathy_lon))
latitude = nc.createVariable('lat', float, ('lat',))
longitude = nc.createVariable('lon', float, ('lon',))
wind_save = nc.createVariable('wind', float, ('lat','lon',))
latitude[:] = bathy_lat
longitude[:] = bathy_lon
wind_save[:] = bathy_data
nc.close()

os.system('gmt grdsample INPUT_bathy.nc -GOUTPUT_bathy.nc -I1.0/1.0 -R-180/180/-79/-60')
os.system('rm INPUT_bathy.nc')

nc = Dataset('OUTPUT_bathy.nc', 'r')
bathy_lat = nc.variables['y'][:]
bathy_lon = nc.variables['x'][:]
bathy_data = nc.variables['z'][:]
nc.close()
os.system('rm OUTPUT_bathy.nc')

# Load the bottom pressure data
M = 1
fig = pl.figure()
for station_name  in  ['DrakePassageSouth_bpr','DrakePassageSouthDeep_bpr','DrakePassageNorthDeep_bpr']:#,'Argentine_Islands_tide_gauge','Syowa_station_tide_gauge']:

    pressure_file = '/Users/jmh2g09/Documents/PhD/Data/BPR/Processed/' + station_name + '_corrected_monthly.csv'

    year = []
    month = []
    dates = []
    bpr = []

    f = open(pressure_file, 'r')
    for line in f:
        line = line.strip()
        columns = line.split(',')
        year.append(int(columns[0]))
        month.append(int(columns[1]))
        # Convert the bottom pressure (millibars) to water depth changes
        # 1 millibar approx = 1 cm water
        if station_name[-3:] == 'bpr':
            bpr.append(float(columns[2]) / 100)
        else:
            bpr.append(float(columns[2]))
        
        dates.append(date(int(columns[0]), int(columns[1]), 15))
        
    f.close()
    
    start_year = year[0]
    start_month = month[0]
    end_year = year[-1]
    end_month = month[-1]
    
    new_year = []
    new_month = []

    for iyear in range(start_year, end_year + 1):
        for imonth in range(1, 13):
            new_year.append(iyear)
            new_month.append(imonth)
    
    new_bpr = np.full(len(new_year), np.NaN)
    for idat in range(len(new_year)):
        if np.any(np.logical_and(np.array(year) == new_year[idat], np.array(month) == new_month[idat])):
            ind = np.where(np.logical_and(np.array(year) == new_year[idat], np.array(month) == new_month[idat]))
            new_bpr[idat] = bpr[ind[0][0]]
    
    ibpr = np.isfinite(new_bpr)
    
    new_bpr = new_bpr[ibpr]
    new_year = np.array(new_year)[ibpr]
    new_month = np.array(new_month)[ibpr]
    dates = np.array(dates)[ibpr]
    
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

    pl.plot(dates, new_bpr - (M * 0.1), label=station_name)
    M += 1
    
    # Load the altimetry data for the timeseries defined by the bottom
    # pressure/tide gauge record

    dot_anom_ts = np.zeros((59, 361, len(new_year)))
    dot_2_anom_ts = np.zeros((59, 361, len(new_year)))

    for t in range(len(new_year)):

        # Make the month string consistent with the file convention
        if int(new_month[t]) < 10:
            month_str = '0' + str(new_month[t])
        elif int(new_month[t]) >= 10:
            month_str = str(new_month[t])

        altimetry_file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + str(new_year[t]) + month_str + '_grid.nc'
        if np.isfinite(new_bpr[t]):
            if os.path.exists(altimetry_file):
                nc = Dataset(altimetry_file, 'r')
                lat = nc.variables['lat'][:]
                lon = nc.variables['lon'][:]
                dot_anom_ts[:, :, t] = nc.variables['dynamic_ocean_topography_anomaly_seasonal_offset'][:]
                dot_2_anom_ts[:, :, t] = nc.variables['dynamic_ocean_topography_anomaly_constant_offset'][:]
                nc.close()
            else:
                dot_anom_ts[:, :, t] = np.full((59, 361), fill_value=np.NaN)
                dot_2_anom_ts[:, :, t] = np.full((59, 361), fill_value=np.NaN)
        else:
            dot_anom_ts[:, :, t] = np.full((59, 361), fill_value=np.NaN)
            dot_2_anom_ts[:, :, t] = np.full((59, 361), fill_value=np.NaN)

    ## Calculate the correlation between the time series and the altimetry data
    dot_anom_xcorr = np.full((59, 361), fill_value=np.NaN)
    dot_anom_xcorr_pvalues = np.full((59, 361), fill_value=np.NaN)
    dot_2_anom_xcorr = np.full((59, 361), fill_value=np.NaN)
    dot_2_anom_xcorr_pvalues = np.full((59, 361), fill_value=np.NaN)

    for ilat in range(len(lat)):
        for ilon in range(len(lon)):
            # If there are nans in the altimetry data
            if np.sum(np.isfinite(dot_anom_ts[ilat, ilon, :])) > len(dot_anom_ts[ilat, ilon, :]) // 2:
                
                its = np.isfinite(dot_anom_ts[ilat, ilon, :])
                
                ts_1 = dot_anom_ts[ilat, ilon, its]
                ts_2 = dot_2_anom_ts[ilat, ilon, its]
                
                new_bpr_2 = new_bpr[its]
                
                xcorr = stats.pearsonr(ts_1, new_bpr_2)
                dot_anom_xcorr[ilat, ilon] = xcorr[0]
                dot_anom_xcorr_pvalues[ilat, ilon] = xcorr[1]
                xcorr_2 = stats.pearsonr(ts_2, new_bpr_2)
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
    
    grid_lats_bathy, grid_lons_bathy = np.meshgrid(bathy_lat, bathy_lon)
    stereo_x_bathy, stereo_y_bathy = m(grid_lons_bathy, grid_lats_bathy)
    
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anom_xcorr)), cmap='RdBu_r')
    m.colorbar()
    pl.clim(1, -1)
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anom_xcorr_pvalues)), [0.05], colors='k')
    
    m.contour(stereo_x_bathy, stereo_y_bathy, np.transpose(np.ma.masked_invalid(-bathy_data)), [3500], colors='g')

    location_x, location_y = m(station_lon, station_lat)
    m.scatter(location_x, location_y, marker='*', s=50, color='k', zorder=100)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/BPR/Figures/' + station_name + '_correlation.png',
        transparent=True, dpi=300, bbox_inches='tight')
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
    m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_2_anom_xcorr_pvalues)), [0.05], colors='k')

    location_x, location_y = m(station_lon, station_lat)
    m.scatter(location_x, location_y, marker='*', s=50, color='k', zorder=100)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/BPR/Figures/' + station_name + '_correlation_constant_offset.png',
        transparent=True, dpi=300, bbox_inches='tight')
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
    
    m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anom_xcorr)) - np.transpose(np.ma.masked_invalid(dot_2_anom_xcorr)), cmap='RdBu_r')
    m.colorbar()
    pl.clim(.1, -.1)
    location_x, location_y = m(station_lon, station_lat)
    m.scatter(location_x, location_y, marker='*', s=50, color='k', zorder=100)
    pl.title('CS-2 and ' + station_name.replace('_', ' ') + ' seasonal - constant offset')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/BPR/Figures/' + station_name + '_correlation_seasonal-constant_offset.png',
        transparent=True, dpi=300, bbox_inches='tight')
    pl.close()

DOT_dates = []
timeseries = []
timeseries_2 = []
f = open('/Users/jmh2g09/Documents/PhD/Data/BPR/Processed/DOT_SAMarea_ts.csv', 'r')
for line in f:
    line = line.strip()
    columns = line.split(' ')
    DOT_dates.append(date(int(columns[0]), int(columns[1]), 15))
    timeseries.append(float(columns[2]))
    timeseries_2.append(float(columns[3]))
f.close()

sam_index = []
f = open('/Users/jmh2g09/Documents/PhD/Data/Wind/SAM_2011-2016.txt', 'r')
for line in f:
    line = line.strip()
    sam_index.append(float(line))
f.close()

pl.plot(DOT_dates, timeseries, label='DOT timeseries (seasonal offset)', ls='--', color='k')
#pl.plot(DOT_dates, np.array(timeseries_2) + 0.1, label='DOT timeseries (no offset)', ls='-.')
#pl.plot(DOT_dates, np.array(timeseries_3) + 0.2, label='DOT timeseries (constant offset)', ls='-.')
pl.plot(DOT_dates, ( - np.array(sam_index) / 100) + 0.1, label='- SAM index / 100')

pl.legend(loc='best', prop={'size':6})
pl.ylabel('Sea level anomaly (m)')
pl.xticks(rotation='vertical')
fig.autofmt_xdate()
#pl.ylim([-0.2, 0.2])
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/BPR/Figures/in_situ_ts.png',
    transparent=True, dpi=300, bbox_inches='tight')
pl.close()