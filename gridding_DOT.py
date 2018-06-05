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

number = np.full((361, 59, 12, 6), fill_value=np.nan)
standard_deviation = np.full((361, 59, 12, 6), fill_value=np.nan)
iciness = np.full((361, 59, 12, 6), fill_value=np.nan)

it = 0
for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        # Cycle through each raw file
        file = '/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + month + '_track.nc'
        print('Gridding: ' + file)
        nc = Dataset(file, 'r')
        lat = nc.variables['latitude'][:]
        lon = nc.variables['longitude'][:]
        dot = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
        dot_2 = nc.variables['dynamic_ocean_topography_constant_offset'][:]
        ssh = nc.variables['sea_surface_height_seasonal_offset'][:]
        ssh_2 = nc.variables['sea_surface_height_constant_offset'][:]
        ice_conc = nc.variables['ice_concentration'][:]
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
        os.system('gmt xyz2grd INPUT_dot.dat -GOUTPUT_n.nc -An -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_dot.dat -GOUTPUT_std.nc -AS -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        
        os.system('gmt xyz2grd INPUT_dot_2.dat -GOUTPUT_dot_2.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_ssh.dat -GOUTPUT_ssh.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_ssh_2.dat -GOUTPUT_ssh_2.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_ice.dat -GOUTPUT_ice.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
        
        os.system('rm INPUT*')
        
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
        
        nc = Dataset('OUTPUT_n.nc', 'r')
        grid_n = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        nc = Dataset('OUTPUT_std.nc', 'r')
        grid_std = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        os.system('rm OUTPUT*')
        
        # Make the longitudes between 0 and 360
        grid_lon[grid_lon < 0] += 361

        grid_dot = grid_dot[np.argsort(grid_lon), :]
        grid_n = grid_n[np.argsort(grid_lon), :]
        grid_std = grid_std[np.argsort(grid_lon), :]
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
            grid_n[land[0][i]][land[1][i]] = np.NaN
            grid_std[land[0][i]][land[1][i]] = np.NaN
            grid_dot_2[land[0][i]][land[1][i]] = np.NaN
            grid_ssh[land[0][i]][land[1][i]] = np.NaN
            grid_ssh_2[land[0][i]][land[1][i]] = np.NaN
        
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
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/Figures/' + year +'/' 
            + year + month + '_DOT.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
        pl.close()
        
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
        m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(grid_n))
        c = m.colorbar()
        c.set_label('Number of Data Points per Node')
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/Figures/' + year +'/' 
            + year + month + '_n.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
        pl.close()
        
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
        m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(grid_std))
        c = m.colorbar()
        c.set_label('Standard Deviation within Node')
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/Figures/' + year +'/' 
            + year + month + '_std.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
        pl.close()
        
        standard_error = grid_std / np.sqrt(grid_n)
        
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
        m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(standard_error))
        c = m.colorbar()
        c.set_label('Standard Deviation within Node')
        pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/Figures/' + year +'/' 
            + year + month + '_error.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
        pl.close()
        
        number[:, :, int(month)-1, int(year)-2011] = grid_n
        standard_deviation[:, :, int(month)-1, int(year)-2011] = grid_std
        iciness[:, :, int(month)-1, int(year)-2011] = grid_ice
        it += 1

        # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Gridded     
        nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc', 'w')

        nc.createDimension('lat', np.size(grid_lat))
        nc.createDimension('lon', np.size(grid_lon))

        latitudes = nc.createVariable('lat', float, ('lat',))
        longitudes = nc.createVariable('lon', float, ('lon',))
        gridded_dot = nc.createVariable('dynamic_ocean_topography_seasonal_offset', float, ('lat','lon'))
        gridded_n = nc.createVariable('number', float, ('lat','lon'))
        gridded_std = nc.createVariable('standard_deviation', float, ('lat','lon'))
        gridded_dot_2 = nc.createVariable('dynamic_ocean_topography_constant_offset', float, ('lat','lon'))
        gridded_ssh = nc.createVariable('sea_surface_height_seasonal_offset', float, ('lat','lon'))
        gridded_ssh_2 = nc.createVariable('sea_surface_height_constant_offset', float, ('lat','lon'))
        gridded_ice = nc.createVariable('ice_concentration', float, ('lat','lon'))

        latitudes.standard_name = 'latitude'
        latitudes.units = 'degrees_north'
        longitudes.standard_name = 'longitude'
        longitudes.units = 'degrees_east'
        gridded_dot.standard_name = 'sea_surface_height_above_EIGEN6c4_seasonal_offset'
        gridded_dot.units = 'm'
        gridded_n.standard_name = 'number_of_data_points_per_node'
        gridded_std.standard_name = 'standard_deviation_in_each_node'
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
        gridded_n[:] = np.transpose(grid_n)
        gridded_std[:] = np.transpose(grid_std)
        gridded_dot_2[:] = np.transpose(grid_dot_2)
        gridded_ssh[:] = np.transpose(grid_ssh)
        gridded_ssh_2[:] = np.transpose(grid_ssh_2)
        gridded_ice[:] = np.transpose(grid_ice)

        nc.close()
        
        print('Complete!')

standard_error = standard_deviation / np.sqrt(number)
ice_standard_error = standard_deviation / np.sqrt(number)
ocean_standard_error = standard_deviation / np.sqrt(number)

ocean_standard_error[iciness >= 1] = np.nan
ice_standard_error[iciness < 1] = np.nan

print('Ice:' + str(np.nanmean(np.nanmean(np.nanmean(np.nanmean(ice_standard_error))))))
print('Ocean:' + str(np.nanmean(np.nanmean(np.nanmean(np.nanmean(ocean_standard_error))))))

# Ice:0.00586079707193
# Ocean:0.00282609141883