import numpy as np
import os
import matplotlib.pyplot as pl
from netCDF4 import Dataset
import functions as funct
from mpl_toolkits.basemap import Basemap
from scipy import stats

# Calculate the average SSH anomaly for every month, resulting in the seasonal cycle 

# Cycle through the years and months
circumpolar_anomaly_cycle = np.zeros(12)
circumpolar_under_ice_cycle = np.zeros(12)
circumpolar_open_ocean_cycle = np.zeros(12)
count = np.zeros(12)

timeseries = []
timeseries_ice = []
timeseries_ocean = []

total_area_ice = []
total_area_ocean = []
total_area_all = []

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/SSH/' + year + '/Anomalies')
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        month_index = int(month) - 1
        # Open the month file 
        # If the file exists...
        file = year + month + '_SSH_anomaly.nc'
        if os.path.isfile(file):
            # Extract the dot anomaly
            nc = Dataset(file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            # Need to load the DOT three times, to avoid calculations crossing over
            ssh_anom  = nc.variables['sea_surface_height_anomaly'][:]
            ssh_anom_ice  = nc.variables['sea_surface_height_anomaly'][:]
            ssh_anom_ocean = nc.variables['sea_surface_height_anomaly'][:]
            ice_conc = nc.variables['sea_ice_concentration'][:]
            nc.close()

            ## Calculate the surface area of each cell
            # Mesh the lat and lon together calculate the surface area for each cell
            grid_lon, grid_lat = np.meshgrid(lon, lat)
            # Calculate the surface area of each cell
            S = funct.surface_area(grid_lat, grid_lon, 0.5, 1.0)
            
            ## Get circumpolar average for all ocean for this month
            # Get total area for the ocean region
            total_area = np.nansum(np.nansum(~np.isnan(ssh_anom) * S))
            total_area_all.append(total_area)
            # Calculate the weighted average DOT for the whole Southern Ocean
            timeseries.append(np.nansum(ssh_anom * S) / total_area)
            circumpolar_anomaly_cycle[month_index] += np.nansum(ssh_anom * S) / total_area

            ## Calculate under ice average for this month
            # Make open ocean and land grid cells NaNs
            open_ocean = np.where(ice_conc == 0)
            for i in range(np.shape(open_ocean)[1]):
                ssh_anom_ice[open_ocean[0][i]][open_ocean[1][i]] = np.nan         
            # Calculate the total area covered by sea ice
            ice_area = np.nansum(np.nansum(~np.isnan(ssh_anom_ice) * S))
            total_area_ice.append(ice_area)
            # Calculate the weighted average DOT under the sea ice
            circumpolar_under_ice_cycle[month_index] += np.nansum(ssh_anom_ice * S) / ice_area
            timeseries_ice.append(np.nansum(ssh_anom_ice * S) / ice_area)
            
            ## Calculate open ocean average for this month
            # Make open ocean and land grid cells NaNs
            under_ice = np.where(ice_conc > 0)
            for i in range(np.shape(under_ice)[1]):
                ssh_anom_ocean[under_ice[0][i]][under_ice[1][i]] = np.nan
            # Calculate the total area covered by sea ice
            ocean_area = np.nansum(np.nansum(~np.isnan(ssh_anom_ocean) * S))
            total_area_ocean.append(ocean_area)
            # Calculate the weighted average DOT under the sea ice
            circumpolar_open_ocean_cycle[month_index] += np.nansum(ssh_anom_ocean * S) / ocean_area
            timeseries_ocean.append(np.nansum(ssh_anom_ocean * S) / ocean_area)

            #test.append((np.nansum(dot_anom_ocean * S) + np.nansum(dot_anom_ice * S)) / total_area)
            count[month_index] += 1
            
            # Cycle through the longitudes, and average the sectors
            #for ilon in np.arange(len(lon)):
             #   if lon[ilon] > 300.:
              #      #print('Weddell')
               #     #weddell_anomaly_cycle[month_index] =+ np.nanmean(dot_anom[ilon, :])
                #    pass
circumpolar_anomaly_cycle /= count
circumpolar_under_ice_cycle /= count
circumpolar_open_ocean_cycle /= count

## Calculate the linear trend for each timeseries
time = np.arange(1, len(timeseries) + 1) / 12

regress_all = stats.linregress(time, timeseries)
line_all = (time * regress_all[0]) + regress_all[1]
print('Total trend: ', regress_all[0]*1000, 'mm / year')

regress_ocean = stats.linregress(time, timeseries_ocean)
line_ocean = (time * regress_ocean[0]) + regress_ocean[1]
print('Ocean trend: ', regress_ocean[0]*1000, 'mm / year')

regress_ice = stats.linregress(time, timeseries_ice)
line_ice = (time * regress_ice[0]) + regress_ice[1]
print('Ice trend: ', regress_ice[0]*1000, 'mm / year')

pl.figure()
pl.plot(total_area_all, label='Total Ocean')
pl.plot(total_area_ice, label='Under Ice')
pl.plot(total_area_ocean, label='Open Ocean')
pl.legend(loc='best')
pl.xlabel('Month')
pl.ylabel('Total surface area (km$^2$)')
pl.savefig('../../../../Seasonal/SSH_total_ice_area.png', format='png', dpi=300)
pl.close()

pl.figure()
pl.plot(np.arange(1, 13), circumpolar_anomaly_cycle, label='Total Ocean')
pl.plot(np.arange(1, 13), circumpolar_under_ice_cycle, label='Under Ice')
pl.plot(np.arange(1, 13), circumpolar_open_ocean_cycle, label='Open Ocean')
pl.legend()
pl.xlabel('Month')
pl.ylabel('Sea Surface Height Anomaly (m)')
pl.xlim([1, 12])
pl.savefig('../../../../Seasonal/SSH_seasonal_cycle.png', format='png', dpi=300)
pl.close()

pl.figure()
pl.plot(timeseries,'k', label='Total Ocean: ' + str((regress_all[0] * 1000).round(2)) + ' mm yr$^{-1}$')
pl.plot(line_all, 'k')
pl.plot(timeseries_ice, 'b', label='Under Ice: ' + str((regress_ice[0] * 1000).round(2)) + ' mm yr$^{-1}$')
pl.plot(line_ice, 'b')
pl.plot(timeseries_ocean, 'r', label='Open Ocean: ' + str((regress_ocean[0] * 1000).round(2)) + ' mm yr$^{-1}$')
pl.plot(line_ocean, 'r')
pl.legend(loc='upper left')
pl.xlabel('Month')
pl.ylabel('Sea Surface Height Anomaly (m)')
pl.savefig('../../../../Seasonal/SSH_seasonal_timeseries.png', format='png', dpi=300)
pl.close()
