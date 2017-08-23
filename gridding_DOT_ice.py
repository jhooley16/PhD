import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import os

lat_resolution = '0.5'    #input('What latitude resolution?: ')
lon_resolution = '1.0'      #input ('What longitude resolution?: ')

os.chdir('/Users/jmh2g09/Documents/PhD/Data/Processed/')

for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        # Cycle through each raw file
        file = '/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + month + '_ice.nc'
        print('Gridding: ' + file)
        nc = Dataset(file, 'r')
        lat = nc.variables['latitude'][:]
        lon = nc.variables['longitude'][:]
        dot = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
        dot_2 = nc.variables['dynamic_ocean_topography_constant_offset'][:]
        ssh = nc.variables['sea_surface_height_seasonal_offset'][:]
        ssh_2 = nc.variables['sea_surface_height_constant_offset'][:]
        ice_conc = nc.variables['sea_ice_concentration'][:]
        nc.close()
        
        input_dot = open('INPUT_dot.dat', 'w')
        input_dot_2 = open('INPUT_dot_2.dat', 'w')
        input_ssh = open('INPUT_ssh.dat', 'w')
        input_ssh_2 = open('INPUT_ssh_2.dat', 'w')
        input_ice = open('INPUT_ice.dat', 'w')
        for ilen in range(len(lat)):
            print(lon[ilen], lat[ilen], dot[ilen], file=input_dot)
            print(lon[ilen], lat[ilen], dot_2[ilen], file=input_dot_2)
            print(lon[ilen], lat[ilen], ssh[ilen], file=input_ssh)
            print(lon[ilen], lat[ilen], ssh_2[ilen], file=input_ssh_2)
            print(lon[ilen], lat[ilen], ice_conc[ilen], file=input_ice)
        input_dot.close()
        input_dot_2.close()
        input_ssh.close()
        input_ssh_2.close()
        input_ice.close()
        
        os.system('gmt xyz2grd INPUT_dot.dat -GOUTPUT_dot.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_dot_2.dat -GOUTPUT_dot_2.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_ssh.dat -GOUTPUT_ssh.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_ssh_2.dat -GOUTPUT_ssh_2.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_ice.dat -GOUTPUT_ice.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        
        os.system('rm INPUT_dot.dat INPUT_dot_2.dat INPUT_ssh.dat INPUT_ssh_2.dat INPUT_ice.dat')
        
        nc = Dataset('OUTPUT_dot.nc', 'r')
        grid_lat = nc.variables['lat'][:]
        grid_lon = nc.variables['lon'][:]
        grid_dot = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        nc = Dataset('OUTPUT_dot_2.nc', 'r')
        grid_dot_2 = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        nc = Dataset('OUTPUT_ssh.nc', 'r')
        grid_ssh = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        nc = Dataset('OUTPUT_ssh_2.nc', 'r')
        grid_ssh_2 = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        nc = Dataset('OUTPUT_ice.nc', 'r')
        grid_ice = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        os.system('rm OUTPUT_dot.nc OUTPUT_dot_2.nc OUTPUT_ssh.nc OUTPUT_ssh_2.nc OUTPUT_ice.nc')

        # Make the longitudes between 0 and 360
        grid_lon[grid_lon < 0] += 361

        grid_dot = grid_dot[np.argsort(grid_lon), :]
        grid_dot_2 = grid_dot_2[np.argsort(grid_lon), :]
        grid_ssh = grid_ssh[np.argsort(grid_lon), :]
        grid_ssh_2 = grid_ssh_2[np.argsort(grid_lon), :]
        grid_ice = grid_ice[np.argsort(grid_lon), :]
        grid_lon = grid_lon[np.argsort(grid_lon)]
        
        # Apply the land mask
        nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/Masks/mask.nc', 'r')
        # Load the mask (ocean == 1)
        ocean_mask = np.transpose(nc.variables['z'][:])
        nc.close()

        land = np.where(ocean_mask != 1)

        for i in range(np.shape(land)[1]):
            grid_dot[land[0][i]][land[1][i]] = np.NaN
            grid_dot_2[land[0][i]][land[1][i]] = np.NaN
            grid_ssh[land[0][i]][land[1][i]] = np.NaN
            grid_ssh_2[land[0][i]][land[1][i]] = np.NaN

        month = file[4:6]
        
        pl.figure()
        pl.clf()
        m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
        m.drawmapboundary()
        m.drawcoastlines(zorder=10)
        m.fillcontinents(zorder=10)
        m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
        grid_lats, grid_lons = np.meshgrid(grid_lat, grid_lon)
        stereo_x, stereo_y = m(grid_lons, grid_lats)
        m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(grid_dot), cmap='RdBu_r')
        c = m.colorbar()
        c.set_label('DOT (m)')
        pl.clim(0, -2.5)
        #pl.clim(np.mean(np.ma.masked_invalid(grid_dot)) - 3*np.std(np.ma.masked_invalid(grid_dot)), np.mean(np.ma.masked_invalid(grid_dot)) + 3*np.std(np.ma.masked_invalid(grid_dot)))
        m.contour(stereo_x, stereo_y, np.ma.masked_invalid(grid_ice), [20,])
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/'+ year +'/Figures/' 
            + year + month + '_DOT_ice.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
        pl.close()

        # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Gridded     
        nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year + month + 'grid_ice.nc', 'w')

        nc.createDimension('lat', np.size(grid_lat))
        nc.createDimension('lon', np.size(grid_lon))

        latitudes = nc.createVariable('latitude', float, ('lat',))
        longitudes = nc.createVariable('longitude', float, ('lon',))
        gridded_dot = nc.createVariable('dynamic_ocean_topography_seasonal_offset', float, ('lat','lon'))
        gridded_dot_2 = nc.createVariable('dynamic_ocean_topography_constant_offset', float, ('lat','lon'))
        gridded_ssh = nc.createVariable('sea_surface_height_seasonal_offset', float, ('lat','lon'))
        gridded_ssh_2 = nc.createVariable('sea_surface_height_constant_offset', float, ('lat','lon'))
        gridded_ice = nc.createVariable('sea_ice_concentration', float, ('lat','lon'))

        latitudes.standard_name = 'latitude'
        latitudes.units = 'degrees_north'
        longitudes.standard_name = 'longitude'
        longitudes.units = 'degrees_east'
        gridded_dot.standard_name = 'sea_surface_height_above_EIGEN6c4_seasonal_offset'
        gridded_dot.units = 'm'
        gridded_dot_2.standard_name = 'sea_surface_height_above_EIGEN6c4_constant_offset'
        gridded_dot_2.units = 'm'
        gridded_ssh.standard_name = 'sea_surface_height_above_WGS84_seasonal_offset'
        gridded_ssh.units = 'm'
        gridded_ssh_2.standard_name = 'sea_surface_height_above_WGS84_constant_offset'
        gridded_ssh_2.units = 'm'
        gridded_ice.standard_name = 'sea_ice_concentration'
        gridded_ice.units = '%'

        latitudes[:] = grid_lat
        longitudes[:] = grid_lon
        gridded_dot[:] = np.transpose(grid_dot)
        gridded_dot_2[:] = np.transpose(grid_dot_2)
        gridded_ssh[:] = np.transpose(grid_ssh)
        gridded_ssh_2[:] = np.transpose(grid_ssh_2)
        gridded_ice[:] = np.transpose(grid_ice)

        nc.close()
        
        print('Complete!')