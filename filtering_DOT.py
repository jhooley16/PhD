import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

for year in ['2011']:#, '2012', '2013', '2014', '2015', '2016']:
    for month in ['01']:#, '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        print(year, month)
        file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc'
        if os.path.isfile(file):
            # Load the lat and lon data for each file
            nc = Dataset(file, 'r')
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            nc.close()

            # Filter the DOT with gaussian filter
            os.system('gmt grdfilter ' + file + '?"dynamic_ocean_topography_seasonal_offset" -D4 -Fg600 -Nr -f0x,1y -GDOT_filt.nc')
            os.system('gmt grdfilter ' + file + '?"dynamic_ocean_topography_constant_offset" -D4 -Fg600 -Nr -f0x,1y -GDOT_2_filt.nc')
            os.system('gmt grdfilter ' + file + '?"sea_surface_height_seasonal_offset" -D4 -Fg600 -Nr -f0x,1y -GSSH_filt.nc')
            os.system('gmt grdfilter ' + file + '?"sea_surface_height_constant_offset" -D4 -Fg600 -Nr -f0x,1y -GSSH_2_filt.nc')
            # Filter the sea ice data
            os.system('gmt grdfilter ' + file + '?"ice_concentration" -D4 -Fg500 -Ni -f0y -f1x -GICE_filt.nc')

            # Open the filtered data
            nc = Dataset('DOT_filt.nc', 'r')
            dot_filt = np.array(nc.variables['z'][:])
            nc.close()
            nc = Dataset('DOT_2_filt.nc', 'r')
            dot_2_filt = np.array(nc.variables['z'][:])
            nc.close()
            nc = Dataset('SSH_filt.nc', 'r')
            ssh_filt = np.array(nc.variables['z'][:])
            nc.close()
            nc = Dataset('SSH_2_filt.nc', 'r')
            ssh_2_filt = np.array(nc.variables['z'][:])
            nc.close()
            nc = Dataset('ICE_filt.nc', 'r')
            ice_filt = np.array(nc.variables['z'][:])
            nc.close()
        
            os.system('rm DOT_filt.nc DOT_2_filt.nc SSH_filt.nc SSH_2_filt.nc ICE_filt.nc')

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
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dot_filt)), cmap='RdBu_r')
            c = m.colorbar()
            pl.clim(0, -2.3)
            m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_filt)), [20,], colors='k')
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/Figures/' + year + '/' + year + month + '_DOT_filt.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
            pl.close()
        
            # Save the gridded data
            os.system('ncks -O -x -v filtered_dynamic_ocean_topography_seasonal_offset,\
filtered_dynamic_ocean_topography_constant_offset,\
filtered_sea_surface_height_seasonal_offset,\
filtered_sea_surface_height_constant_offset,\
filtered_sea_ice_concentration ' + file + ' ' + file)
            nc = Dataset(file, 'a')
            dot = nc.createVariable('filtered_dynamic_ocean_topography_seasonal_offset', float, ('lat','lon'))
            dot_2 = nc.createVariable('filtered_dynamic_ocean_topography_constant_offset', float, ('lat','lon'))
            ssh = nc.createVariable('filtered_sea_surface_height_seasonal_offset', float, ('lat','lon'))
            ssh_2 = nc.createVariable('filtered_sea_surface_height_constant_offset', float, ('lat','lon'))
            sea_ice_concentration = nc.createVariable('filtered_sea_ice_concentration', float, ('lat', 'lon'))
        
            dot.standard_name = 'filtered_sea_surface_height_above_EIGEN6c4_seasonal_offset'
            dot.units = 'm'
            dot_2.standard_name = 'filtered_sea_surface_height_above_EIGEN6c4_constant_offset'
            dot_2.units = 'm'
            ssh.standard_name = 'filtered_sea_surface_height_above_WGS84_seasonal_offset'
            ssh.units = 'm'
            ssh_2.standard_name = 'filtered_sea_surface_height_above_WGS84_constant_offset'
            ssh_2.units = 'm'
            sea_ice_concentration.standard_name = 'filtered_sea_ice_concentration'
            sea_ice_concentration.units = '%'
        
            dot[:] = dot_filt
            dot_2[:] = dot_2_filt
            ssh[:] = ssh_filt
            ssh_2[:] = ssh_2_filt
            sea_ice_concentration[:] = ice_filt
            nc.close()