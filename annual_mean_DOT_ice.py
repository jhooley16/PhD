import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

# Initiate the arrays to be filled, separating out the years and the months
dot_array = np.full((59, 361, 12, 7), fill_value=np.NaN)
dot_array_2 = np.full((59, 361, 12, 7), fill_value=np.NaN)
ssh_array = np.full((59, 361, 12, 7), fill_value=np.NaN)
ssh_array_2 = np.full((59, 361, 12, 7), fill_value=np.NaN)
ice_array = np.full((59, 361, 12, 7), fill_value=np.NaN)

# Cycle through the months and arrange the grids by year and month
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid_ice.nc'
        if os.path.isfile(file):
            nc = Dataset(file, 'r')
            dot_array[:, :, int(month)-1, int(year)-2010] = nc.variables['filtered_dynamic_ocean_topography_seasonal_offset'][:]
            dot_array_2[:, :, int(month)-1, int(year)-2010] = nc.variables['filtered_dynamic_ocean_topography_constant_offset'][:]
            ssh_array[:, :, int(month)-1, int(year)-2010] = nc.variables['filtered_sea_serface_height_seasonal_offset'][:]
            ssh_array_2[:, :, int(month)-1, int(year)-2010] = nc.variables['filtered_sea_serface_height_constant_offset'][:]
            ice_array[:, :, int(month)-1, int(year)-2010] = nc.variables['filtered_sea_ice_concentration'][:]
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            nc.close()

# Take an average of all years and all months for a total mean
dot_mean = np.nanmean(np.nanmean(dot_array, axis=3), axis=2)
dot_2_mean = np.nanmean(np.nanmean(dot_array_2, axis=3), axis=2)
ssh_mean = np.nanmean(np.nanmean(ssh_array, axis=3), axis=2)
ssh_2_mean = np.nanmean(np.nanmean(ssh_array_2, axis=3), axis=2)

# Take the mean over only the years, leaving a monthly mean for each
dot_mean_months = np.nanmean(dot_array, axis=3)
dot_mean_months_2 = np.nanmean(dot_array_2, axis=3)
ssh_mean_months = np.nanmean(ssh_array, axis=3)
ssh_mean_months_2 = np.nanmean(ssh_array_2, axis=3)
# Do the same for the ice, producing a grid of mean ice concentration
ice_mean_months = np.nanmean(ice_array, axis=3)

# Save the all month mean
nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/Means/mean_ice.nc', 'w')

nc.createDimension('lat', np.size(grid_lat))
nc.createDimension('lon', np.size(grid_lon))
nc.createDimension('months', 12)

latitudes = nc.creatVariable('latitude', float, ('lat',))
longitudes = nc.creatVariable('longitude', float, ('lon',))
dot_annual_mean = nc.createVariable('mean_dynamic_topography_seasonal_offset', float, ('lat','lon'))
dot_2_annual_mean = nc.createVariable('mean_dynamic_topography_constant_offset', float, ('lat','lon'))
ssh_annual_mean = nc.createVariable('mean_sea_surface_height_seasonal_offset', float, ('lat','lon'))
ssh_2_annual_mean = nc.createVariable('mean_sea_surface_height_constant_offset', float, ('lat','lon'))
dot_monthly_mean = nc.createVariable('monthly_mean_dynamic_topography_seasonal_offset', float, ('lat','lon', 'months'))
dot_2_monthly_mean = nc.createVariable('monthly_mean_dynamic_topography_constant_offset', float, ('lat','lon', 'months'))
ssh_monthly_mean = nc.createVariable('monthly_sea_surface_height_constant_offset', float, ('lat','lon', 'months'))
ssh_2_monthly_mean = nc.createVariable('monthly_mean_sea_surface_height_constant_offset', float, ('lat','lon', 'months'))
ice_monthly_mean = nc.createVariable('monthly_mean_sea_ice_concentration', float, ('lat','lon', 'months'))

dot_annual_mean.standard_name = 'mean_sea_surface_height_above_EIGEN6c4_seasonal_offset'
dot_annual_mean.units = 'm'
dot_2_annual_mean.standard_name = 'mean_sea_surface_height_above_EIGEN6c4_constant_offset'
dot_2_annual_mean.units = 'm'
ssh_annual_mean.standard_name = 'mean_sea_surface_height_above_WGS84_seasonal_offset'
ssh_annual_mean.units = 'm'
ssh_2_annual_mean.standard_name = 'mean_sea_surface_height_above_WGS84_constant_offset'
ssh_2_annual_mean.units = 'm'
dot_monthly_mean.standard_name = 'monthly_mean_sea_surface_height_above_EIGEN6c4_seasonal_offset'
dot_monthly_mean.units = 'm'
dot_2_monthly_mean.standard_name = 'monthly_mean_sea_surface_height_above_EIGEN6c4_constant_offset'
dot_2_monthly_mean.units = 'm'
ssh_monthly_mean.standard_name = 'monthly_mean_sea_surface_height_above_WGS84_seasonal_offset'
ssh_monthly_mean.units = 'm'
ssh_2_monthly_mean.standard_name = 'monthly_mean_sea_surface_height_above_WGS84_constant_offset'
ssh_2_monthly_mean.units = 'm'
ice_monthly_mean.standard_name = 'monthly_mean_sea_ice_concentration'
ice_monthly_mean.units = '%'

dot_annual_mean[:] = dot_mean
dot_2_annual_mean[:] = dot_2_mean
ssh_annual_mean[:] = ssh_mean
ssh_2_annual_mean[:] = ssh_2_mean
dot_monthly_mean[:] = dot_mean_months
dot_2_monthly_mean[:] = dot_mean_months_2
ssh_monthly_mean[:] = ssh_mean_months
ssh_2_monthly_mean[:] = ssh_mean_months_2
ice_monthly_mean[:] = ice_mean_months
nc.close()