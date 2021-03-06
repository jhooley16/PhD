import os
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import shiftgrid

## take a zonal average of the DOT around the continent 

# Load the mean DOT calculated earlier...
nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/Means/mean.nc', 'r')
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
dot  = nc.variables['mean_dynamic_topography_seasonal_offset'][:]
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
pl.ylim([-2.2, -2])
pl.xlim([0, 1000])
pl.grid()
#pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('MDT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/mean_zonal_slope.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
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
    for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
        # Change directory and see if the file exists...
        file_ssh = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc'
        if os.path.isfile(file_ssh):
            nc = Dataset(file_ssh, 'r')
            lat_season = nc.variables['lat'][:]
            lon_season = nc.variables['lon'][:]
            dot_season = nc.variables['filtered_dynamic_ocean_topography_seasonal_offset'][:]
            u_vel_season = nc.variables['surface_u_velocity'][:]
            nc.close()
            
            # Find longest string of data in each meridional section
            length_of_values_season = []
            for ilon in lon_season[:270]:
                ilon = int(ilon)
                length_of_values_season.append(np.nansum(np.isfinite(dot_season[:, ilon])))

            new_grid_season = np.full((max(length_of_values_season), 360), fill_value=np.NaN)
            new_u_grid_season = np.full((max(length_of_values_season), 360), fill_value=np.NaN)  
            # pad out difference with NaNs
            for ilon in lon_season[:270]:
                ilon = int(ilon)
                difference_season_dot = max(length_of_values_season) - len(dot_season[np.isfinite(dot_season[:, ilon]), ilon])
                difference_season_u = max(length_of_values_season) - len(u_vel_season[np.isfinite(u_vel_season[:, ilon]), ilon])
                new_grid_season[:, ilon] = np.append(dot_season[np.isfinite(dot_season[:, ilon]), ilon], [np.NaN] * difference_season_dot)
                new_u_grid_season[:, ilon] = np.append(u_vel_season[np.isfinite(u_vel_season[:, ilon]), ilon], [np.NaN] * difference_season_u)
                
            # Take the average in the zonal direction
            zonal_mean_ssh_season[:len(np.nanmean(new_grid_season, axis=1)), it] = np.nanmean(new_grid_season, axis=1)
            zonal_mean_u_season[:len(np.nanmean(new_u_grid_season, axis=1)), it] = np.nanmean(new_u_grid_season, axis=1)
#         Change directory and see if the file exists...
#         os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/velocity/' + year)
#         file_vel = year + month + '_velocity.nc'
#         if os.path.isfile(file_vel):
#             nc = Dataset(file_vel, 'r')
#             lat_season = nc.variables['latitude'][:]
#             lon_season_u = nc.variables['longitude'][:]
#             u_vel_season = nc.variables['surface_u_velocity'][:]
#             nc.close()
#             
#             u_vel_season, lon_season_u = shiftgrid(300.,pre_u_vel_season,pre_lon_season,start=False)
#             
#             Find longest string of data in each meridional section
#             length_of_values_season = []
#             for ilon in lon_season_u[:270]:
#                 ilon = int(ilon)
#                 length_of_values_season.append(np.nansum(np.isfinite(u_vel_season[:, ilon])))
#             
#             new_u_grid_season = np.full((max(length_of_values_season), 360), fill_value=np.NaN)
#             pad out difference with NaNs
#             for ilon in lon_season_u[:270]:
#                 ilon = int(ilon)
#                 difference_season = max(length_of_values_season) - len(u_vel_season[np.isfinite(u_vel_season[:, ilon]), ilon])
#                 new_u_grid_season[:, ilon] = np.append(u_vel_season[np.isfinite(u_vel_season[:, ilon]), ilon], [np.NaN] * difference_season)
#             Take the average in the zonal direction
#             zonal_mean_u_season[:len(np.nanmean(new_u_grid_season, axis=1)), it] = np.nanmean(new_u_grid_season, axis=1)

        it += 1
    # Take the average for all the months
    monthly_zonal_ssh_average[:, it_mean] = np.nanmean(zonal_mean_ssh_season, axis=1)
    monthly_zonal_u_average[:, it_mean] = np.nanmean(zonal_mean_u_season, axis=1)
    it_mean += 1

months = ['Jan', 'Feb', 'Mar','Apr', 'May', 'Jun', 'July', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']

faster_zonal_ssh_average = np.nanmean(monthly_zonal_ssh_average[:, 4:7], axis=1)
slower_zonal_ssh_average_1 = np.nanmean(monthly_zonal_ssh_average[:, :2], axis=1)
slower_zonal_ssh_average_2 = np.nanmean(monthly_zonal_ssh_average[:, 11:], axis=1)
slower_zonal_ssh_average = (slower_zonal_ssh_average_1 + slower_zonal_ssh_average_2) / 2

pl.figure()
pl.plot(distance, faster_zonal_ssh_average, label='AMJ', linestyle='-', color='k')
pl.scatter(distance, faster_zonal_ssh_average, color='k')
pl.plot(distance, slower_zonal_ssh_average, label='NDJF', linestyle='--', color='k')
pl.scatter(distance, slower_zonal_ssh_average, color='k')
pl.legend(loc='best')
pl.grid()
pl.xlim([0, 600])
pl.ylim([-2.2, -2])
#pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('MDT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/zonal_slope_compare.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

faster_zonal_u_average = np.nanmean(monthly_zonal_u_average[:, 2:7], axis=1)
slower_zonal_u_average_1 = np.nanmean(monthly_zonal_u_average[:, :2], axis=1)
slower_zonal_u_average_2 = np.nanmean(monthly_zonal_u_average[:, 7:], axis=1)
slower_zonal_u_average = (slower_zonal_u_average_1 + slower_zonal_u_average_2) / 2

pl.figure()
pl.plot(distance, faster_zonal_u_average, label='AMJ', linestyle='-', color='k')
pl.scatter(distance, faster_zonal_u_average, color='k')
pl.plot(distance, slower_zonal_u_average, label='NDJF', linestyle='--', color='k')
pl.scatter(distance, slower_zonal_u_average, color='k')
pl.legend(loc='best')
pl.grid()
pl.xlim([0, 600])
pl.ylim([-0.06, 0.06])
#pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/zonal_u_compare.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
for T in range(12):
    pl.plot(distance, monthly_zonal_ssh_average[:, T], label=months[T])
pl.legend(loc='lower right')
pl.xlim([0, 600])
pl.ylim([-2.2, -2])
#pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('MDT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/full_zonal_slope.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
for T in range(12):
    pl.plot(distance, monthly_zonal_u_average[:, T], label=months[T])
pl.plot(distance, [0] * len(distance), linestyle='--', color='k')
pl.legend(loc='best')
pl.xlim([0, 600])
pl.ylim([-0.06, 0.06])
#pl.title('Zonal u-velocity mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/full_zonal_u.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
for T in range(2, 7):
    pl.plot(distance, monthly_zonal_ssh_average[:, T], label=months[T])
pl.legend(loc='best')
pl.grid()
pl.xlim([0, 600])
pl.ylim([-2.2, -1.9])
pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('MDT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/winter_zonal_slope.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
for T in range(2, 7):
    pl.plot(distance, monthly_zonal_u_average[:, T], label=months[T])
pl.legend(loc='best')
pl.grid()
pl.xlim([0, 600])
pl.ylim([-0.02, 0.03])
pl.title('u-velocity mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/winter_zonal_u.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
for T in range(7, 12):
    pl.plot(distance, monthly_zonal_ssh_average[:, T], label=months[T])
for T in range(2):
    pl.plot(distance, monthly_zonal_ssh_average[:, T], label=months[T])
pl.legend(loc='best')
pl.grid()
pl.xlim([0, 600])
pl.ylim([-2.2, -1.9])
pl.title('Zonal DOT mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('MDT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/summer_zonal_slope.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
for T in range(7, 12):
    pl.plot(distance, monthly_zonal_u_average[:, T], label=months[T])
for T in range(2):
    pl.plot(distance, monthly_zonal_u_average[:, T], label=months[T])
pl.legend(loc='best')
pl.grid()
pl.xlim([0, 600])
pl.ylim([-0.02, 0.03])
pl.title('u-velocity mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/summer_zonal_u.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

# Take the gradient of the continental slope slope for each monthly average
grad = np.full(12, fill_value=np.NaN)
for month in range(12):
    grad[month] = np.nanmean(np.gradient(monthly_zonal_ssh_average[:, month] * 1000, np.mean(np.gradient(distance)))[:4])

pl.figure()
pl.plot(range(1, 13), grad, linestyle='-', marker='o', color='k')
pl.xlim([0.5, 12.5])
pl.grid()
pl.xticks(np.arange(1, 13), months)
pl.ylabel('Meridional Gradient (mm / km )')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/monthly_zonal_slope_gradient.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

# Take the average u-velocity
u_velocity_mean = np.nanmean(monthly_zonal_u_average, axis=1)

pl.figure()
pl.plot(distance, u_velocity_mean, linestyle='-', marker='o')
pl.grid()
pl.xlim([0, 600])
pl.ylim([-0.06, 0.06])
#pl.title('Zonal u-velocity mean between 0$^\circ$E - 270$^\circ$E')
pl.xlabel('Distance from Antarctic Coast (km)')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/mean_zonal_u.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

# Take the gradient of the continental slope slope for each monthly average
grad_u = np.full(12, fill_value=np.NaN)
for month in range(12):
    grad_u[month] = np.nanmean(monthly_zonal_u_average[:, month])

pl.figure()
pl.plot(range(1, 13), grad_u, linestyle='-', marker='o')
pl.xlim([0.5, 12.5])
pl.title('Average u-velocity of slope along Antarctic Coastline')
pl.xlabel('Month')
pl.ylabel('u-velocity (m s$^{-1}$)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/monthly_zonal_u.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()