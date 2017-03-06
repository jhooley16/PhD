import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import os

for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    lat_resolution = '0.5'    #input('What latitude resolution?: ')
    lon_resolution = '1.0'      #input ('What longitude resolution?: ')

    # Look at the raw data files
    os.chdir('/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + '/DOT_track')

    # Cycle through each raw file
    for file in os.listdir():
        if file[-12:] == 'DOT_track.nc':
            print('Gridding: ' + file)
            nc = Dataset(file, 'r')
            lat = nc.variables['latitude'][:]
            lon = nc.variables['longitude'][:]
            dot = nc.variables['dynamic_ocean_topography_seasonal_offset'][:]
            dot_2 = nc.variables['dynamic_ocean_topography_no_offset'][:]
            dot_3 = nc.variables['dynamic_ocean_topography_constant_offset'][:]
            ice_conc = nc.variables['sea_ice_concentration'][:]
            nc.close()
            
            input_dot = open('INPUT_dot.dat', 'w')
            input_dot_2 = open('INPUT_dot_2.dat', 'w')
            input_dot_3 = open('INPUT_dot_3.dat', 'w')
            input_ice = open('INPUT_ice.dat', 'w')
            for ilen in range(len(lat)):
                print(lon[ilen], lat[ilen], dot[ilen], file=input_dot)
                print(lon[ilen], lat[ilen], dot_2[ilen], file=input_dot_2)
                print(lon[ilen], lat[ilen], dot_3[ilen], file=input_dot_3)
                print(lon[ilen], lat[ilen], ice_conc[ilen], file=input_ice)
            input_dot.close()
            input_ice.close()
            input_dot_2.close()
            input_dot_3.close()
            
            os.system('gmt xyz2grd INPUT_dot.dat -GOUTPUT_dot.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
            os.system('gmt xyz2grd INPUT_dot_2.dat -GOUTPUT_dot_2.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
            os.system('gmt xyz2grd INPUT_dot_3.dat -GOUTPUT_dot_3.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
            os.system('gmt xyz2grd INPUT_ice.dat -GOUTPUT_ice.nc -I' + lon_resolution + '/' + lat_resolution + ' -R-180/180/-79/-50 -fig')
            
            os.system('rm INPUT_dot.dat INPUT_dot_2.dat INPUT_dot_3.dat INPUT_ice.dat')
            
            nc = Dataset('OUTPUT_dot.nc', 'r')
            grid_lat = nc.variables['lat'][:]
            grid_lon = nc.variables['lon'][:]
            grid_dot = np.array(np.transpose(nc.variables['z'][:]))
            nc.close()
            
            nc = Dataset('OUTPUT_dot_2.nc', 'r')
            grid_dot_2 = np.array(np.transpose(nc.variables['z'][:]))
            nc.close()
            
            nc = Dataset('OUTPUT_dot_3.nc', 'r')
            grid_dot_3 = np.array(np.transpose(nc.variables['z'][:]))
            nc.close()
            
            nc = Dataset('OUTPUT_ice.nc', 'r')
            grid_ice = np.array(np.transpose(nc.variables['z'][:]))
            nc.close()
            
            os.system('rm OUTPUT_dot.nc OUTPUT_dot_2.nc OUTPUT_dot_3.nc OUTPUT_ice.nc')

            # Make the longitudes between 0 and 360
            grid_lon[grid_lon < 0] += 361

            grid_dot = grid_dot[np.argsort(grid_lon), :]
            grid_dot_2 = grid_dot_2[np.argsort(grid_lon), :]
            grid_dot_3 = grid_dot_3[np.argsort(grid_lon), :]
            grid_ice = grid_ice[np.argsort(grid_lon), :]
            grid_lon = grid_lon[np.argsort(grid_lon)]
            
            # Apply the land mask
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/mask.nc', 'r')
            # Load the mask (ocean == 1)
            ocean_mask = np.transpose(nc.variables['z'][:])
            nc.close()

            land = np.where(ocean_mask != 1)

            for i in range(np.shape(land)[1]):
                grid_dot[land[0][i]][land[1][i]] = np.NaN
                grid_dot_2[land[0][i]][land[1][i]] = np.NaN
                grid_dot_3[land[0][i]][land[1][i]] = np.NaN

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
                + year + month + '_DOT_gridded_.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
            pl.close()
    
            # Put the data in a .nc file in /Users/jmh2g09/Documents/PhD/Data/Gridded     
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/Gridded/DOT/' + year + '/' + year + month + '_DOT.nc', 'w')

            nc.createDimension('lat', np.size(grid_lat))
            nc.createDimension('lon', np.size(grid_lon))

            latitudes = nc.createVariable('latitude', float, ('lat',))
            longitudes = nc.createVariable('longitude', float, ('lon',))
            gridded_dot = nc.createVariable('dynamic_ocean_topography_seasonal_offset', float, ('lat','lon'))
            gridded_dot_2 = nc.createVariable('dynamic_ocean_topography_no_offset', float, ('lat','lon'))
            gridded_dot_3 = nc.createVariable('dynamic_ocean_topography_constant_offset', float, ('lat','lon'))
            gridded_ice = nc.createVariable('sea_ice_concentration', float, ('lat','lon'))

            latitudes.long_name = 'latitude'
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.long_name = 'longitude'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            gridded_dot.long_name = 'dynamic_ocean_topography'
            gridded_dot.standard_name = 'sea_surface_height_above_EIGEN6c4_geoid_seasonal_retracker_offset_correction'
            gridded_dot.units = 'm'
            gridded_dot_2.long_name = 'dynamic_ocean_topography'
            gridded_dot_2.standard_name = 'sea_surface_height_above_EIGEN6c4_geoid_no_retracker_offset_correction'
            gridded_dot_2.units = 'm'
            gridded_dot_3.long_name = 'dynamic_ocean_topography'
            gridded_dot_3.standard_name = 'sea_surface_height_above_EIGEN6c4_geoid_constant_retracker_offset_correction'
            gridded_dot_3.units = 'm'
            gridded_ice.long_name = 'sea_ice_concentration'
            gridded_ice.standard_name = 'sea_ice_concentration'
            gridded_ice.units = '%'

            latitudes[:] = grid_lat
            longitudes[:] = grid_lon
            gridded_dot[:] = np.transpose(grid_dot)
            gridded_dot_2[:] = np.transpose(grid_dot_2)
            gridded_dot_3[:] = np.transpose(grid_dot_3)
            gridded_ice[:] = np.transpose(grid_ice)

            nc.close()
            
            if year == '2011':
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
                m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(grid_dot) - np.ma.masked_invalid(grid_dot_2), cmap='RdBu_r')
                m.colorbar()
                pl.title('Seasonal Offset: ' + month)
                #pl.clim(0, -2.5)
                #pl.clim(np.mean(np.ma.masked_invalid(grid_dot)) - 3*np.std(np.ma.masked_invalid(grid_dot)), np.mean(np.ma.masked_invalid(grid_dot)) + 3*np.std(np.ma.masked_invalid(grid_dot)))
                m.contour(stereo_x, stereo_y, np.ma.masked_invalid(grid_ice), [20,])
                pl.savefig('/Users/jmh2g09/Documents/PhD/Data/SeparateModes/Figures/offset_seasonal_' + month + '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
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
                m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(grid_dot_3) - np.ma.masked_invalid(grid_dot_2), cmap='RdBu_r')
                m.colorbar()
                pl.title('Constant Offset: ' + month)
                #pl.clim(0, -2.5)
                #pl.clim(np.mean(np.ma.masked_invalid(grid_dot)) - 3*np.std(np.ma.masked_invalid(grid_dot)), np.mean(np.ma.masked_invalid(grid_dot)) + 3*np.std(np.ma.masked_invalid(grid_dot)))
                m.contour(stereo_x, stereo_y, np.ma.masked_invalid(grid_ice), [20,])
                pl.savefig('/Users/jmh2g09/Documents/PhD/Data/SeparateModes/Figures/offset_constant_' + month + '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
                pl.close()
            
            print('Complete!')