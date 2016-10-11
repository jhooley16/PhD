import os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import shiftgrid

## take a zonal average of the DOT around the continent 

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/')

# Load the mean DOT calculated earlier...
nc = Dataset('MDT_mean.nc', 'r')
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
dot  = nc.variables['mean_dynamic_topography'][:]
nc.close()

## Need to rearrange data so that grid cells next to the AA coast are first
# Find the longest length of a meridional section that has data (58 in this case)
length_of_values = []
for ilon in lon[:271]:
    length_of_values.append(np.sum(np.isfinite(dot[:, ilon])))

# Make a new grid with a new meridional length
new_grid = np.full((60, 270), fill_value=np.NaN)

# For every meridional section, pad out the difference with nans if the data 
# isn't long enough
for ilon in lon[:270]:
    ilon = int(ilon)
    difference = max(length_of_values) - len(dot[np.isfinite(dot[:, ilon]), ilon])
    new_grid[:len(np.append(dot[np.isfinite(dot[:, ilon]), ilon], [np.NaN] * difference)), ilon] = np.append(dot[np.isfinite(dot[:, ilon]), ilon], [np.NaN] * difference)

# Take the mean in the zonal (x) direction
zonal_mean = np.nanmean(new_grid, axis=1)
# Get the distance between points (latitude, so all equal), 
# 0.5 increment cells away from the coast
distance = np.arange(0.5, 60 // 2 + 0.5,  0.5) * 60 * 1.862 # km

pl.figure()
pl.plot(distance, zonal_mean, linestyle='-', marker='o')
pl.ylim([-1.8, -1.6])
pl.xlim([0, 1000])
pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/mean_zonal_slope.png', format='png', transparent=True, dpi=300)
pl.close()

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
            
            #dot_season, lon_season = shiftgrid(300.,pre_dot_season,pre_lon_season,start=False)
            
            
            # Find longest string of data in each meridional section
            length_of_values_season = []
            for ilon in lon_season[:270]:
                ilon = int(ilon)
                length_of_values_season.append(np.nansum(np.isfinite(dot_season[:, ilon])))

            new_grid_season = np.full((max(length_of_values_season), 360), fill_value=np.NaN)
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
            lon_season_u = nc.variables['longitude'][:]
            u_vel_season = nc.variables['surface_u_velocity'][:]
            nc.close()
            
            #u_vel_season, lon_season_u = shiftgrid(300.,pre_u_vel_season,pre_lon_season,start=False)
            
            # Find longest string of data in each meridional section
            length_of_values_season = []
            for ilon in lon_season_u[:270]:
                ilon = int(ilon)
                length_of_values_season.append(np.nansum(np.isfinite(u_vel_season[:, ilon])))
            
            new_u_grid_season = np.full((max(length_of_values_season), 360), fill_value=np.NaN)
            # pad out difference with NaNs
            for ilon in lon_season_u[:270]:
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

months = ['Jan', 'Feb', 'Mar','Apr', 'May', 'Jun', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']

pl.figure()
for T in range(12):
    pl.plot(distance, monthly_zonal_ssh_average[:, T], label=months[T])
pl.legend(loc='lower right')
pl.xlim([0, 1000])
pl.ylim([-2.5, -2])
pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/full_zonal_slope.png', format='png', transparent=True, dpi=300)
pl.close()

pl.figure()
for T in range(12):
    pl.plot(distance, monthly_zonal_u_average[:, T], label=months[T])
pl.legend(loc='lower right')
pl.xlim([0, 1000])
pl.ylim([-2.5, -2])
pl.title('Zonal u-velocity mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/full_zonal_u.png', format='png', transparent=True, dpi=300)
pl.close()

pl.figure()
for T in range(4):
    pl.plot(distance, monthly_zonal_ssh_average[:, T], label=months[T])
pl.legend()
pl.xlim([0, 1000])
pl.ylim([-2.5, -2])
pl.title('Summer Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/summer_zonal_slope.png', format='png', transparent=True, dpi=300)
pl.close()

pl.figure()
for T in range(4):
    pl.plot(distance, monthly_zonal_u_average[:, T], label=months[T])
pl.legend()
pl.xlim([0, 1000])
pl.ylim([-2.5, -2])
pl.title('Summer u-velocity mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/summer_zonal_u.png', format='png', transparent=True, dpi=300)
pl.close()

pl.figure()
for T in range(5, 12):
    pl.plot(distance, monthly_zonal_ssh_average[:, T], label=months[T])
pl.legend()
pl.xlim([0, 1000])
#pl.ylim([-1.8, -1.6])
pl.title('Winter Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/winter_zonal_slope.png', format='png', transparent=True, dpi=300)
pl.close()

pl.figure()
for T in range(5, 12):
    pl.plot(distance, monthly_zonal_u_average[:, T], label=months[T])
pl.legend()
pl.xlim([0, 1000])
#pl.ylim([-1.8, -1.6])
pl.title('Winter u-velocity mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/winter_zonal_u.png', format='png', transparent=True, dpi=300)
pl.close()

# Take the gradient of the continental slope slope for each monthly average
grad = np.full(12, fill_value=np.NaN)
for month in range(12):
    grad[month] = np.nanmean(np.gradient(monthly_zonal_ssh_average[:, month] * 1000, np.mean(np.gradient(distance)))[:4])

pl.figure()
pl.plot(range(1, 13), grad, linestyle='-', marker='o')
pl.xlim([0.5, 12.5])
pl.title('Average gradient of slope along Antarctic Coastline')
pl.xlabel('Month')
pl.ylabel('Sea Surface Gradient (mm / km )')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/monthly_zonal_slope_gradient.png', format='png', transparent=True, dpi=300)
pl.close()

# Take the average u-velocity
u_velocity_mean = np.nanmean(monthly_zonal_u_average, axis=1)

pl.figure()
pl.plot(distance, u_velocity_mean, linestyle='-', marker='o')
pl.xlim([0, 1000])
pl.ylim([-2.5, -2])
pl.title('Zonal u-velocity mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/mean_zonal_u.png', format='png', transparent=True, dpi=300)
pl.close()

# Take the gradient of the continental slope slope for each monthly average
grad_u = np.full(12, fill_value=np.NaN)
for month in range(12):
    grad_u[month] = np.nanmean(np.gradient(monthly_zonal_u_average[:, month], np.mean(np.gradient(distance)))[:4])

pl.figure()
pl.plot(range(1, 13), grad, linestyle='-', marker='o')
pl.xlim([0.5, 12.5])
pl.title('Average u-velocity of slope along Antarctic Coastline')
pl.xlabel('Month')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Coastal Current/monthly_zonal_u.png', format='png', transparent=True, dpi=300)
pl.close()