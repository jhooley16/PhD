import os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import shiftgrid

## Find the zonal average for every month
# initiate the array and the counter
monthly_zonal_ssh_average = np.full((60, 12), fill_value=np.NaN)
monthly_zonal_u_average = np.full((60, 12), fill_value=np.NaN)
it_mean = 0
# Cycle through each month and take the average for all the Jans, Febs, etc...
for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
    zonal_mean_ssh_season = np.full((60, 7), fill_value=np.NaN)
    zonal_mean_u_season = np.full((60, 7), fill_value=np.NaN)
    it = 0
    # For each year, take the average for all Jans, Febs, etc...
    for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
        # Change directory and see if the file exists...
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
        file_ssh = year + month + '_DOT.nc'
        if os.path.isfile(file_ssh):
            nc = Dataset(file_ssh, 'r')
            lat_season = nc.variables['latitude'][:]
            lon_season = nc.variables['longitude'][:]
            dot_season = nc.variables['dynamic_ocean_topography'][:]
            nc.close()
            
            dotgrid, newlon = shiftgrid(300.,dot_season,lon_season,start=False)
            
            # Find longest string of data in each meridional section
            length_of_values_season = []
            for ilon in lon_season[:271]:
                ilon = int(ilon)
                length_of_values_season.append(np.nansum(np.isfinite(dot_season[:, ilon])))
            
            new_grid_season = np.full((max(length_of_values_season), 270), fill_value=np.NaN)
            # pad out difference with NaNs
            for ilon in lon_season[:270]:
                ilon = int(ilon)
                difference_season = max(length_of_values_season) - len(dot_season[np.isfinite(dot_season[:, ilon]), ilon])
                new_grid_season[:, ilon] = np.append(dot_season[np.isfinite(dot_season[:, ilon]), ilon], [np.NaN] * difference_season)
                
            # Take the average in the zonal direction
            zonal_mean_ssh_season[:len(np.nanmean(new_grid_season, axis=1)), it] = np.nanmean(new_grid_season, axis=1)

        # Change directory and see if the file exists...
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/velocity/' + year)
        file_vel = year + month + '_velocity.nc'
        if os.path.isfile(file_vel):
            nc = Dataset(file_vel, 'r')
            lat_season = nc.variables['latitude'][:]
            lon_season = nc.variables['longitude'][:]
            u_vel_season = nc.variables['surface_u_velocity'][:]
            nc.close()
            # Find longest string of data in each meridional section
            length_of_values_season = []
            for ilon in lon_season[:271]:
                ilon = int(ilon)
                length_of_values_season.append(np.nansum(np.isfinite(u_vel_season[:, ilon])))
            
            new_u_grid_season = np.full((max(length_of_values_season), 270), fill_value=np.NaN)
            # pad out difference with NaNs
            for ilon in lon_season[:270]:
                ilon = int(ilon)
                difference_season = max(length_of_values_season) - len(u_vel_season[np.isfinite(u_vel_season[:, ilon]), ilon])
                new_u_grid_season[:, ilon] = np.append(u_vel_season[np.isfinite(u_vel_season[:, ilon]), ilon], [np.NaN] * difference_season)
            # Take the average in the zonal direction
            zonal_mean_u_season[:len(np.nanmean(new_u_grid_season, axis=1)), it] = np.nanmean(new_u_grid_season, axis=1)

        it += 1
    # Take the average for all the months
    monthly_zonal_ssh_average[:, it_mean] = np.nanmean(zonal_mean_ssh_season, axis=1)
    monthly_zonal_u_average[:, it_mean] = np.nanmean(zonal_mean_u_season, axis=1)
    it_mean += 1