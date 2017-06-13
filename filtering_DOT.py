import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year)
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if os.path.isfile(year + month + '_DOT.nc'):
            print(year, month)
            # Load the lat and lon data for each file
            nc = Dataset(year + month + '_DOT.nc', 'r')
    
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]

            nc.close()

            # Filter the DOT with a 500km gaussian filter
            os.system('gmt grdfilter ' + year + month + '_DOT.nc?"dynamic_ocean_topography_seasonal_offset" -D4 -Fg1000 -Nr -f0y -f1x -GDOT_filt.nc')
            os.system('gmt grdfilter ' + year + month + '_DOT.nc?"dynamic_ocean_topography_no_offset" -D4 -Fg1000 -Nr -f0y -f1x -GDOT_2_filt.nc')
            os.system('gmt grdfilter ' + year + month + '_DOT.nc?"dynamic_ocean_topography_constant_offset" -D4 -Fg1000 -Nr -f0y -f1x -GDOT_3_filt.nc')
            
            # Open the filtered data
            nc = Dataset('DOT_filt.nc', 'r')
            dot_filt = np.array(nc.variables['z'][:])
            nc.close()
            
            nc = Dataset('DOT_2_filt.nc', 'r')
            dot_2_filt = np.array(nc.variables['z'][:])
            nc.close()
            
            nc = Dataset('DOT_3_filt.nc', 'r')
            dot_3_filt = np.array(nc.variables['z'][:])
            nc.close()

            os.system('gmt grdfilter ' + year + month + '_DOT.nc?"sea_ice_concentration" -D4 -Fg500 -Ni -f0y -f1x -GICE_filt.nc')

            nc = Dataset('ICE_filt.nc', 'r')
            ice_filt = np.array(nc.variables['z'][:])
            nc.close()

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
            #pl.clim(np.mean(np.ma.masked_invalid(dot_filt)) + 3*np.std(np.ma.masked_invalid(dot_filt)), np.mean(np.ma.masked_invalid(dot_filt)) - 3*np.std(np.ma.masked_invalid(dot_filt)))
            m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_filt)), [20,], colors='k')
            pl.savefig('Figures/' + year + month + '_DOT_filt.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
            pl.close()

            nc = Dataset(year + month + '_DOT_filt.nc', 'w')

            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))

            latitudes = nc.createVariable('latitude', float, ('lat',))
            longitudes = nc.createVariable('longitude', float, ('lon',))
            sea_surface_height = nc.createVariable('dynamic_ocean_topography_seasonal_offset', float, ('lat','lon'))
            sea_surface_height_2 = nc.createVariable('dynamic_ocean_topography_no_offset', float, ('lat','lon'))
            sea_surface_height_3 = nc.createVariable('dynamic_ocean_topography_constant_offset', float, ('lat','lon'))
            sea_ice_concentration = nc.createVariable('sea_ice_concentration', float, ('lat', 'lon'))
        
            latitudes[:] = lat
            longitudes[:] = lon
            sea_surface_height[:] = dot_filt
            sea_surface_height_2[:] = dot_2_filt
            sea_surface_height_3[:] = dot_3_filt
            sea_ice_concentration[:] = ice_filt
            nc.close()
            
            os.system('rm DOT_filt.nc DOT_2_filt.nc DOT_3_filt.nc ICE_filt.nc')