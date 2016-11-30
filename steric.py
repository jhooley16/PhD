import numpy as np
from netCDF4 import Dataset
from datetime import date
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import functions as funct

steric_ts = []
dot_ts = []
eustatic_ts = []

steric_ts_ice = []
dot_ts_ice = []
eustatic_ts_ice = []

steric_ts_ocean = []
dot_ts_ocean = []
eustatic_ts_ocean = []

dates = []

total_area_ice = []
total_area_ice = []

steric_seasonal = np.full((12, 5), fill_value=np.NaN)
dot_seasonal = np.full((12, 5), fill_value=np.NaN)
eustatic_seasonal =np.full((12, 5), fill_value=np.NaN)

steric_seasonal_ice = np.full((12, 5), fill_value=np.NaN)
dot_seasonal_ice = np.full((12, 5), fill_value=np.NaN)
eustatic_seasonal_ice =np.full((12, 5), fill_value=np.NaN)

steric_seasonal_ocean = np.full((12, 5), fill_value=np.NaN)
dot_seasonal_ocean = np.full((12, 5), fill_value=np.NaN)
eustatic_seasonal_ocean =np.full((12, 5), fill_value=np.NaN)

trigger = False
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    
    # Cycle through the months
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if year == '2010' and month == '11':
            trigger = True
        if year == '2016' and month == '03':
            trigger = False
        
        if trigger == True:
            dates.append(date(int(year), int(month), 15))
            
            DOT_file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year + '/Anomalies/' + year + month + '_DOT_anomaly.nc'
            # Load DOT
            nc = Dataset(DOT_file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            dot = nc.variables['dynamic_ocean_topography_anomaly_seasonal_offset'][:] # (lat, lon)
            dot_ice = nc.variables['dynamic_ocean_topography_anomaly_seasonal_offset'][:] # (lat, lon)
            dot_ocean = nc.variables['dynamic_ocean_topography_anomaly_seasonal_offset'][:] # (lat, lon)
            ice = nc.variables['sea_ice_concentration'][:]
            nc.close()
            
            GRACE_file = '/Users/jmh2g09/Documents/PhD/Data/GRACE/' + year + month + '_GRACE.nc'
            # Load GRACE
            nc = Dataset(GRACE_file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            eustatic = nc.variables['GRACE'][:] # (lat, lon)
            eustatic_ice = nc.variables['GRACE'][:] # (lat, lon)
            eustatic_ocean = nc.variables['GRACE'][:] # (lat, lon)
            nc.close()
            
            # Calculate the steric 
            steric = dot - eustatic
            steric_ice = dot_ice - eustatic_ice
            steric_ocean = dot_ocean - eustatic_ocean

            ## Calculate the surface area of each cell
            # Mesh the lat and lon together calculate the surface area for each cell
            grid_lon, grid_lat = np.meshgrid(lon, lat)
            # Calculate the surface area of each cell
            S = funct.surface_area(grid_lat, grid_lon, 0.5, 1.0)
            # Get total area for the ocean region
            total_area = np.nansum(np.nansum(~np.isnan(dot) * S))
            
            # Define the open ocean as 0% ice concentration
            open_ocean = np.where(ice == 0)
            under_ice = np.where(ice > 0)
            
            for i in range(np.shape(open_ocean)[1]):
                steric_ice[open_ocean[0][i]][open_ocean[1][i]] = np.NaN
                dot_ice[open_ocean[0][i]][open_ocean[1][i]] = np.NaN
                eustatic_ice[open_ocean[0][i]][open_ocean[1][i]] = np.NaN
            
            for i in range(np.shape(under_ice)[1]):
                steric_ocean[under_ice[0][i]][under_ice[1][i]] = np.NaN
                dot_ocean[under_ice[0][i]][under_ice[1][i]] = np.NaN
                eustatic_ocean[under_ice[0][i]][under_ice[1][i]] = np.NaN
                
            # Calculate the total area covered by sea ice
            ice_area = np.nansum(np.nansum(~np.isnan(dot_ice) * S))
            # Calculate the total area NOT covered by sea ice
            ocean_area = np.nansum(np.nansum(~np.isnan(dot_ocean) * S))
            
            # Calculate weighed mean for all the ocean
            steric_ts.append(np.nansum(steric * S) / total_area)
            dot_ts.append(np.nansum(dot * S) / total_area)
            eustatic_ts.append(np.nansum(eustatic * S) / total_area)
            
            # Calculate the weighted mean for the ice-covered ocean
            steric_ts_ice.append(np.nansum(steric_ice * S) / ice_area)
            dot_ts_ice.append(np.nansum(dot_ice * S) / ice_area)
            eustatic_ts_ice.append(np.nansum(eustatic_ice * S) / ice_area)
            
            # Calculate the weighted mean for the open ocean
            steric_ts_ocean.append(np.nansum(steric_ocean * S) / ocean_area)
            dot_ts_ocean.append(np.nansum(dot_ocean * S) / ocean_area)
            eustatic_ts_ocean.append(np.nansum(eustatic_ocean * S) / ocean_area)
            
            if int(year) >= 2011 and int(year) < 2016:
                it = int(year) - 2011
                # Calculate the monthly contribution for all the ocean
                steric_seasonal[int(month) - 1, it] = np.nansum(steric * S) / total_area
                dot_seasonal[int(month) - 1, it] = np.nansum(dot * S) / total_area
                eustatic_seasonal[int(month) - 1, it] = np.nansum(eustatic * S) / total_area
                
                # Calculate the monthly contribution for open ocean
                steric_seasonal_ocean[int(month) - 1, it] = np.nansum(steric_ocean * S) / ocean_area
                dot_seasonal_ocean[int(month) - 1, it] = np.nansum(dot_ocean * S) / ocean_area
                eustatic_seasonal_ocean[int(month) - 1, it] = np.nansum(eustatic_ocean * S) / ocean_area
                
                # Calculate the monthly contribution for under ice
                steric_seasonal_ice[int(month) - 1, it] = np.nansum(steric_ice * S) / ice_area
                dot_seasonal_ice[int(month) - 1, it] = np.nansum(dot_ice * S) / ice_area
                eustatic_seasonal_ice[int(month) - 1, it] = np.nansum(eustatic_ice * S) / ice_area


fig = pl.figure()
pl.plot(dates, dot_ts, label='DOT')
pl.plot(dates, np.array(eustatic_ts) + .1, label='Eustatic DOT')
pl.plot(dates, np.array(steric_ts) - .1, label='Steric DOT')
pl.legend(loc='best')
pl.xticks(rotation='vertical')
fig.autofmt_xdate()
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Steric/steric_timeseries.png', transparent=True, doi=300)
pl.close()

fig = pl.figure()
pl.plot(dates, dot_ts_ocean, label='DOT')
pl.plot(dates, np.array(eustatic_ts_ocean) + .1, label='Eustatic DOT')
pl.plot(dates, np.array(steric_ts_ocean) - .1, label='Steric DOT')
pl.legend(loc='best')
pl.xticks(rotation='vertical')
fig.autofmt_xdate()
pl.ylabel('DOT (m)')
pl.title('Open Ocean')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Steric/steric_timeseries_ocean.png', transparent=True, doi=300)
pl.close()

fig = pl.figure()
pl.plot(dates, dot_ts_ice, label='DOT')
pl.plot(dates, np.array(eustatic_ts_ice) + .1, label='Eustatic DOT')
pl.plot(dates, np.array(steric_ts_ice) - .1, label='Steric DOT')
pl.legend(loc='best')
pl.xticks(rotation='vertical')
fig.autofmt_xdate()
pl.ylabel('DOT (m)')
pl.title('Under Ice')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Steric/steric_timeseries_ice.png', transparent=True, doi=300)
pl.close()

steric_season = np.nanmean(steric_seasonal, axis=1)
dot_season = np.nanmean(dot_seasonal, axis=1)
eustatic_season = np.nanmean(eustatic_seasonal, axis=1)

steric_season_ice = np.nanmean(steric_seasonal_ice, axis=1)
dot_season_ice = np.nanmean(dot_seasonal_ice, axis=1)
eustatic_season_ice = np.nanmean(eustatic_seasonal_ice, axis=1)

steric_season_ocean = np.nanmean(steric_seasonal_ocean, axis=1)
dot_season_ocean = np.nanmean(dot_seasonal_ocean, axis=1)
eustatic_season_ocean = np.nanmean(eustatic_seasonal_ocean, axis=1)

fig = pl.figure()
pl.plot(range(1, 13), np.array(dot_season), label='DOT')
pl.plot(range(1, 13), np.array(eustatic_season), label='Eustatic DOT')
pl.plot(range(1, 13), np.array(steric_season), label='Steric DOT')
pl.xlim([1, 12])
pl.legend(loc='best')
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Steric/steric_seasonal.png', transparent=True, doi=300)
pl.close()

fig = pl.figure()
pl.plot(range(1, 13), np.array(dot_season_ice), label='DOT')
pl.plot(range(1, 13), np.array(eustatic_season_ice), label='Eustatic DOT')
pl.plot(range(1, 13), np.array(steric_season_ice), label='Steric DOT')
pl.xlim([1, 12])
pl.legend(loc='best')
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Steric/steric_seasonal_ice.png', transparent=True, doi=300)
pl.close()

fig = pl.figure()
pl.plot(range(1, 13), np.array(dot_season_ocean), label='DOT')
pl.plot(range(1, 13), np.array(eustatic_season_ocean), label='Eustatic DOT')
pl.plot(range(1, 13), np.array(steric_season_ocean), label='Steric DOT')
pl.xlim([1, 12])
pl.legend(loc='best')
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Steric/steric_seasonal_ocean.png', transparent=True, doi=300)
pl.close()