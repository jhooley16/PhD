import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import functions as funct
from scipy.signal import coherence
from datetime import date

# Load the SSHA from the ice-only region and calculate it's coherence with the
# whole Southern Ocean data

ssha1_cat = np.full((64, 59, 361), fill_value=np.NaN)
ssha2_cat = np.full((64, 59, 361), fill_value=np.NaN)
ssha3_cat = np.full((64, 59, 361), fill_value=np.NaN)
ssha_ice_ts = np.full(64, fill_value=np.NaN)
ssha_ice_ts_2 = np.full(64, fill_value=np.NaN)
ssha_ice_ts_3 = np.full(64, fill_value=np.NaN)
DOT_dates = []
cat = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        path = '/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year + '/Anomalies/'
        
        file = year + month + '_DOT_anomaly.nc'
        filepath = path + file
        
        file_ice = year + month + '_DOT_anomaly_ice.nc'
        filepath_ice = path + file_ice
        
        if os.path.isfile(filepath):
        
            nc = Dataset(filepath, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            ssha1_cat[cat, :, :] = nc.variables['dynamic_ocean_topography_anomaly_seasonal_offset'][:]
            ssha2_cat[cat, :, :] = nc.variables['dynamic_ocean_topography_anomaly_constant_offset'][:]
            ssha3_cat[cat, :, :] = nc.variables['dynamic_ocean_topography_anomaly_no_offset'][:]
            nc.close()
            
            nc = Dataset(filepath_ice, 'r')
            ssha_ice_1 = nc.variables['dynamic_ocean_topography_anomaly_seasonal_offset'][:]
            ssha_ice_2 = nc.variables['dynamic_ocean_topography_anomaly_no_offset'][:]
            ssha_ice_3 = nc.variables['dynamic_ocean_topography_anomaly_constant_offset'][:]
            nc.close()
            
            DOT_dates.append(date(int(year), int(month), 15))
            
            ## Calculate the surface area of each cell
            # Mesh the lat and lon together calculate the surface area for each cell
            grid_lon, grid_lat = np.meshgrid(lon, lat)
            # Calculate the surface area of each cell
            S = funct.surface_area(grid_lat, grid_lon, 0.5, 1.0)
            total_area_ice = np.nansum(S * ~np.isnan(ssha_ice_1))
            
            ssha_ice_ts[cat] = np.nansum(ssha_ice_1 * S) / total_area_ice
            ssha_ice_ts_2[cat] = np.nansum(ssha_ice_2 * S) / total_area_ice
            ssha_ice_ts_3[cat] = np.nansum(ssha_ice_3 * S) / total_area_ice
            cat += 1

coherr = np.full((59, 361, 7), fill_value=np.NaN)
coherr_2 = np.full((59, 361, 7), fill_value=np.NaN)
coherr_3 = np.full((59, 361, 7), fill_value=np.NaN)
for ilat in range(len(lat)):
    for ilon in range(len(lon)):
        if np.any(np.isfinite(ssha3_cat[:, ilat, ilon])):
            A = ssha1_cat[:, ilat, ilon]
            B = ssha2_cat[:, ilat, ilon]
            C = ssha3_cat[:, ilat, ilon]
            
            f, Cxy = coherence(ssha_ice_ts, funct.inpaint_nans(A), nperseg=12)
            coherr[ilat, ilon, :] = Cxy
            
            f, Cxy = coherence(ssha_ice_ts_2, funct.inpaint_nans(A), nperseg=12)
            coherr_2[ilat, ilon, :] = Cxy
            
            f, Cxy = coherence(ssha_ice_ts_3, funct.inpaint_nans(A), nperseg=12)
            coherr_3[ilat, ilon, :] = Cxy

station_lon = [320, 30, 135, 190]
station_lat = [10, 22, 28, 5]

fig = pl.figure()
pl.plot(DOT_dates, ssha_ice_ts, label='permanant ice-region average (seasonal)', color='b', ls='-')
pl.plot(DOT_dates, ssha_ice_ts_3, label='permanant ice-region average (constant)', color='b', ls='--')

pl.plot(DOT_dates, ssha1_cat[:, station_lat[0], station_lon[0]] - 0.1, label='1 (-0.1)', color='k', ls='-')
pl.plot(DOT_dates, ssha3_cat[:, station_lat[0], station_lon[0]] - 0.1, label=None, color='k', ls='--')

pl.plot(DOT_dates, ssha1_cat[:, station_lat[1], station_lon[1]] - 0.2, label='2 (-0.2)', color='r', ls='-')
pl.plot(DOT_dates, ssha3_cat[:, station_lat[1], station_lon[1]] - 0.2, label=None, color='r', ls='--')

pl.plot(DOT_dates, ssha1_cat[:, station_lat[2], station_lon[2]] + 0.1, label='3 (+0.1)', color='g', ls='-')
pl.plot(DOT_dates, ssha3_cat[:, station_lat[2], station_lon[2]] + 0.1, label=None, color='g', ls='--')

pl.plot(DOT_dates, ssha1_cat[:, station_lat[3], station_lon[3]] + 0.2, label='4 (+0.2)', color='y', ls='-')
pl.plot(DOT_dates, ssha3_cat[:, station_lat[3], station_lon[3]] + 0.2, label=None, color='y', ls='--')

pl.xticks(rotation='vertical')
fig.autofmt_xdate()
pl.legend(loc='best', prop={'size':6})
pl.ylabel('SSH anomaly (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coherence/timeseries.png', dpi=300, transparent=True)
pl.close()

for num in range(7): 
    pl.figure()
    pl.clf()
    m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    stereo_x, stereo_y = m(grid_lon, grid_lat)
    m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(coherr[:, :, num]), cmap='RdBu_r')
    m.colorbar()
    pl.clim([0, 1])
    location_x, location_y = m(lon[station_lon], lat[station_lat])
    m.scatter(location_x[0], location_y[0], marker='$1$', s=60, color='k', zorder=100)
    m.scatter(location_x[1], location_y[1], marker='$2$', s=60, color='k', zorder=100)
    m.scatter(location_x[2], location_y[2], marker='$3$', s=60, color='k', zorder=100)
    m.scatter(location_x[3], location_y[3], marker='$4$', s=60, color='k', zorder=100)
    pl.title('Coherence at period: ' + str(np.round(1/f[num], 2)) + ' months')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coherence/coherence_' + str(np.round(1/f[num], 2)) + '_months_seasonal_offset.png', dpi=300, transparent=True)
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    stereo_x, stereo_y = m(grid_lon, grid_lat)
    m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(coherr_2[:, :, num]), cmap='RdBu_r')
    m.colorbar()
    pl.clim([0, 1])
    m.scatter(location_x[0], location_y[0], marker='$1$', s=60, color='k', zorder=100)
    m.scatter(location_x[1], location_y[1], marker='$2$', s=60, color='k', zorder=100)
    m.scatter(location_x[2], location_y[2], marker='$3$', s=60, color='k', zorder=100)
    m.scatter(location_x[3], location_y[3], marker='$4$', s=60, color='k', zorder=100)
    pl.title('Coherence at period: ' + str(np.round(1/f[num], 2)) + ' months')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coherence/coherence_' + str(np.round(1/f[num], 2)) + '_months_no_offset.png')
    pl.close()

    pl.figure()
    pl.clf()
    m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    stereo_x, stereo_y = m(grid_lon, grid_lat)
    m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(coherr_3[:, :, num]), cmap='RdBu_r')
    m.colorbar()
    pl.clim([0, 1])
    m.scatter(location_x[0], location_y[0], marker='$1$', s=60, color='k', zorder=100)
    m.scatter(location_x[1], location_y[1], marker='$2$', s=60, color='k', zorder=100)
    m.scatter(location_x[2], location_y[2], marker='$3$', s=60, color='k', zorder=100)
    m.scatter(location_x[3], location_y[3], marker='$4$', s=60, color='k', zorder=100)
    pl.title('Coherence at period: ' + str(np.round(1/f[num], 2)) + ' months')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coherence/coherence_' + str(np.round(1/f[num], 2)) + '_months_constant_offset.png')
    pl.close()