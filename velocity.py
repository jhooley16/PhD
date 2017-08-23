import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
import functions as funct

## For each month, calculate the zonal geostrophic velocity
# Define the constants for the geostrophic equation
dy = 0.5 * 60 * 1862 # distance between two latitude points (0.5 grid spacing) (m)
g = 9.81 # gravitational acceleration (meters per square second)
omega = 2 * np.pi / (24 * 60 * 60) # rotation rate of earth (radians per second)

# define the constants for the longitude distance equation
a = 6378137. # radius of earth (m)

test_u = np.full((59, 361, 72), fill_value=np.NaN)
test_v = np.full((59, 361, 72), fill_value=np.NaN)
test_dot = np.full((59, 361, 72), fill_value=np.NaN)
it = 0
for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc'
        # If the file exists, load the dot data
        print(year, month)
        nc = Dataset(file, 'r')
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        dot  = nc.variables['dynamic_ocean_topography_seasonal_offset'][:] # meters
        nc.close()
        
        # Calculate the longitudinal distance
        dx = (np.pi / 180) * a * np.cos(lat * np.pi / 180)
        dx_2d = np.transpose(np.tile(dx, (361, 1)))

        # Calculate the coriolis parameter
        f = 2 * omega * np.sin(lat * np.pi / 180)
        f_2d = np.transpose(np.tile(f, (361, 1)))

        # Calculate the dot gradients in the x and y directions
        dssh_dy = np.gradient(dot, axis=0) / dy
        
        dssh_dx = np.gradient(dot, axis=1) / dx_2d

        # Calculate the geostrophic currents
        u =  - (g / f_2d) * dssh_dy
        v = (g / f_2d) * dssh_dx
        
        test_u[:, :, it] = - (g / f_2d) * dssh_dy
        test_v[:, :, it] = (g / f_2d) * dssh_dx
        test_dot[:, :, it] = dot
        it += 1
        # Save the data to a file
        os.system('ncks -O -x -v surface_u_velocity,\
surface_v_velocity ' + file + ' ' + file)
        nc = Dataset(file, 'a')

        u_velocity = nc.createVariable('surface_u_velocity', float, ('lat','lon'))
        v_velocity = nc.createVariable('surface_v_velocity', float, ('lat','lon'))

        u_velocity.standard_name = 'surface_zonal_velocity_x'
        u_velocity.units = 'meters_per_second'
        v_velocity.standard_name = 'surface_zonal_velocity_y'
        v_velocity.units = 'meters_per_second'

        u_velocity[:] = u
        v_velocity[:] = v

        nc.close()
            
#             pl.figure()
#             pl.clf()
#             m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
#             m.drawmapboundary()
#             m.drawcoastlines(zorder=10)
#             m.fillcontinents(zorder=10)
#             m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
#             m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
#         
#             grid_lats, grid_lons = np.meshgrid(lat, lon)
#             stereo_x, stereo_y = m(grid_lons, grid_lats)
#         
#             m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(v)), cmap='RdBu_r')
#             m.colorbar()
#             pl.clim(.5, -.5)
#             #pl.clim(np.mean(np.ma.masked_invalid(grid_dot)) - 3*np.std(np.ma.masked_invalid(grid_dot)), np.mean(np.ma.masked_invalid(grid_dot)) + 3*np.std(np.ma.masked_invalid(grid_dot)))
#             m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_conc)), [40,])
#             pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/velocity/'+ year +'/Figures/' 
#                 + year + month + '_v_velocity_.png', format='png', transparent=True, dpi=300)
#             pl.close()
#             
#             pl.figure()
#             pl.clf()
#             m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
#             m.drawmapboundary()
#             m.drawcoastlines(zorder=10)
#             m.fillcontinents(zorder=10)
#             m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
#             m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
#         
#             grid_lats, grid_lons = np.meshgrid(lat, lon)
#             stereo_x, stereo_y = m(grid_lons, grid_lats)
#         
#             m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(u)), cmap='RdBu_r')
#             m.colorbar()
#             pl.clim(.5, -.5)
#             #pl.clim(np.mean(np.ma.masked_invalid(grid_dot)) - 3*np.std(np.ma.masked_invalid(grid_dot)), np.mean(np.ma.masked_invalid(grid_dot)) + 3*np.std(np.ma.masked_invalid(grid_dot)))
#             m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(ice_conc)), [40,])
#             pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/velocity/'+ year +'/Figures/' 
#                 + year + month + '_u_velocity_.png', format='png', transparent=True, dpi=300)
#             pl.close()
#             
#             pl.figure()
#             pl.clf()
#             m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='h')
#             m.drawmapboundary()
#             m.drawcoastlines(zorder=10)
#             m.fillcontinents(zorder=10)
#             m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
#             m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
#         
#             ugrid, newlons = shiftgrid(180.,u,lon,start=False)
#             vgrid, newlons = shiftgrid(180.,v,lon,start=False)
#             dotgrid, newlons = shiftgrid(180.,dot,lon,start=False)
#             icegrid, newlons = shiftgrid(180.,ice_conc,lon,start=False)
#             
#             uproj,vproj,xx,yy = m.transform_vector(ugrid,vgrid,newlons,lat,70, 70,returnxy=True,masked=True)
#             
#             grid_lats, grid_lons = np.meshgrid(lat, newlons)
#             stereo_x, stereo_y = m(grid_lons, grid_lats)
#             
#             m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(dotgrid)), cmap='RdBu_r')
#             m.colorbar()
#             pl.clim(0, -2)
#             m.quiver(xx, yy,  -uproj, -vproj)
#             
#             m.contour(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(icegrid)),  [40,])
#             pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Gridded/velocity/'+ year +'/Figures/' 
#                 + year + month + '_velocity_vectors.png', format='png', transparent=True, dpi=300)
#             pl.close()