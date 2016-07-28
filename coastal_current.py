import os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl

## take a zonal average of the DOT around the continent 

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/')

nc = Dataset('MDT_mean.nc', 'r')
            
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]

dot  = nc.variables['mean_dynamic_topography'][:]
nc.close()

## 
length_of_values = []
for ilon in lon[:271]:
    length_of_values.append(np.sum(np.isfinite(dot[:, ilon])))

new_grid = np.full((max(length_of_values), 270), fill_value=np.NaN)

for ilon in lon[:270]:
    difference = max(length_of_values) - len(dot[np.isfinite(dot[:, ilon]), ilon])
    
    new_grid[:, ilon] = np.append(dot[np.isfinite(dot[:, ilon]), ilon], [np.NaN] * difference)

zonal_mean = np.nanmean(new_grid, axis=1)
distance = np.arange(0.5, max(length_of_values) + 0.5) * 60 * 1.862

pl.figure()
pl.plot(distance, zonal_mean)
pl.xlim([0, 6000])
pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('DOT (m)')
#pl.show()
pl.close()

monthly_zonal_average = np.full((58, 12), fill_value=np.NaN)
it_mean = 0
for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
    zonal_mean_season = np.full((58, 7), fill_value=np.NaN)
    it = 0
    for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
        
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
        file = year + month + '_DOT.nc'
        if os.path.isfile(file):
            # Extract the dot
            nc = Dataset(file, 'r')
            lat_season = nc.variables['latitude'][:]
            lon_season = nc.variables['longitude'][:]

            dot_season = nc.variables['dynamic_ocean_topography'][:]
            nc.close()
            
            length_of_values_season = []
            for ilon in lon_season[:271]:
                length_of_values_season.append(np.sum(np.isfinite(dot_season[:, ilon])))
            
            new_grid_season = np.full((max(length_of_values_season), 270), fill_value=np.NaN)
            
            for ilon in lon_season[:270]:
                difference_season = max(length_of_values_season) - len(dot_season[np.isfinite(dot_season[:, ilon]), ilon])
    
                new_grid_season[:, ilon] = np.append(dot[np.isfinite(dot_season[:, ilon]), ilon], [np.NaN] * difference_season)

            zonal_mean_season[:, it] = np.nanmean(new_grid_season, axis=1)

        it += 1
    monthly_zonal_average[:, it_mean] = np.nanmean(zonal_mean_season, axis=1)
    it_mean += 1

print(np.shape(monthly_zonal_average))

pl.figure()
for T in range(4):
    pl.plot(distance, monthly_zonal_average[:, T])
pl.xlim([0, 2000])
pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('DOT (m)')
pl.show()
pl.close()