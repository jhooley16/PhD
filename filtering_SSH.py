import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/SSH/' + year)
    
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:

        if os.path.isfile(year + month + '_SSH.nc'):
            print(year, month)
            # Load the lat and lon data for each file
            nc = Dataset(year + month + '_SSH.nc', 'r')
    
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]

            nc.close()

            # Filter the DOT with a 300km gaussian filter
            os.system('gmt grdfilter ' + year + month + '_SSH.nc?"sea_surface_height" -D4 -Fg500 -Nr -f0y -f1x -G' + year + month + '_SSH_filt.nc')

            # Open the filtered data
            nc = Dataset(year + month + '_SSH_filt.nc', 'r')
            ssh_filt = np.array(nc.variables['z'][:])
            nc.close()

            os.system('gmt grdfilter ' + year + month + '_SSH.nc?"sea_ice_concentration" -D4 -Fg600 -Ni -f0y -f1x -G' + year + month + '_SSH_filt.nc')

            nc = Dataset(year + month + '_SSH_filt.nc', 'r')
            ice_filt = np.array(nc.variables['z'][:])
            nc.close()

            # Apply the land mask

            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/mask.nc', 'r')
            # Load the mask (ocean == 1)
            ocean_mask = nc.variables['z'][:]
            nc.close()
            
            land = np.where(ocean_mask != 1)
            for i in range(np.shape(land)[1]):
                ssh_filt[land[0][i]][land[1][i]] = np.NaN

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
        
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ssh_filt)), cmap='RdBu_r')
            m.colorbar()
            #pl.clim([0, -2])
            #pl.clim(np.mean(np.ma.masked_invalid(dot_filt)) + 3*np.std(np.ma.masked_invalid(dot_filt)), np.mean(np.ma.masked_invalid(dot_filt)) - 3*np.std(np.ma.masked_invalid(dot_filt)))
            m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_filt)), [60,])
            pl.savefig('Figures/' + year + month + '_SSH_filt.png', format='png', transparent=True, dpi=300)
            pl.close()

            nc = Dataset(year + month + '_SSH_filt.nc', 'w')

            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))

            latitudes = nc.createVariable('latitude', float, ('lat',))
            longitudes = nc.createVariable('longitude', float, ('lon',))
            sea_surface_height = nc.createVariable('sea_surface_height', float, ('lat','lon'))
            sea_ice_concentration = nc.createVariable('sea_ice_concentration', float, ('lat', 'lon'))
        
            latitudes[:] = lat
            longitudes[:] = lon
            sea_surface_height[:] = ssh_filt
            sea_ice_concentration[:] = ice_filt
            nc.close()