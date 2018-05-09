import numpy as np
import os
import matplotlib.pyplot as pl
from netCDF4 import Dataset
import functions as funct
from mpl_toolkits.basemap import Basemap

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

os.system('gmt grdsample INPUT_bathy.nc -GOUTPUT_bathy.nc -I1.0/1.0 -R-180/180/-79/-50')
os.system('rm INPUT_bathy.nc')

nc = Dataset('OUTPUT_bathy.nc', 'r')
bathy_lat = nc.variables['y'][:]
bathy_lon = nc.variables['x'][:]
bathy_data = nc.variables['z'][:]
nc.close()
os.system('rm OUTPUT_bathy.nc')


dot_anom_AMJ = np.full((59, 361, 72), fill_value=np.NaN)
dot_anom_NDJF = np.full((59, 361, 72), fill_value=np.NaN)
t = 0
# Cycle through the years
for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        # If the file exists, open the file and extract the data
        file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc'
        if os.path.isfile(file):
            nc = Dataset(file, 'r')
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            # Assign the monthly gridded DOT anomaly to the array
            if int(month) == 4 or int(month) == 5 or int(month) == 6 or int(month) == 7 or int(month) == 8:
                dot_anom_AMJ[:, :, t]  = nc.variables['dynamic_ocean_topography_anomaly_seasonal_offset'][:]
            elif int(month) == 1 or int(month) == 2 or int(month) == 11 or int(month) == 12:
                dot_anom_NDJF[:, :, t]  = nc.variables['dynamic_ocean_topography_anomaly_seasonal_offset'][:]
            nc.close()
            t += 1

dot_anomaly = np.nanmean(dot_anom_AMJ, 2) - np.nanmean(dot_anom_NDJF, 2)

# Load mean sea ice data

ice_file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/Means/mean.nc'

nc = Dataset(ice_file, 'r')
ice_mean = nc.variables['monthly_mean_sea_ice_concentration'][:]
nc.close()

min_ice = ice_mean[:, :, 2]
max_ice = ice_mean[:, :, 8]

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

bathy_grid_lats, bathy_grid_lons = np.meshgrid(bathy_lat, bathy_lon)
bathy_stereo_x, bathy_stereo_y = m(bathy_grid_lons, bathy_grid_lats)
        
m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_anomaly)), cmap='RdBu_r')
m.colorbar()
pl.clim(0.1, -0.1)
m.contour(bathy_stereo_x, bathy_stereo_y, np.transpose(-bathy_data), [2000, ], colors='g')
m.contour(stereo_x, stereo_y, np.transpose(max_ice), [20, ], colors='k')
m.contour(stereo_x, stereo_y, np.transpose(-min_ice), [-20, ], colors='k')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seasonal/DOT_AMJNDJF_anomaly.png', transparent=True, dpi=300)
pl.close()