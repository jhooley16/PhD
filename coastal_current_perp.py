import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Polygon, Point, LineString
from geopy.distance import vincenty
import functions as funct

m = Basemap(projection='spstere', boundinglat=-60, lon_0=180, resolution='l')

length = 200
dot_array = np.full((6, 12, length, 11), fill_value=np.NaN)
dist_2 = np.full((6, 12, length, 11), fill_value=np.NaN)
coast_distance_array = np.full((6, 12, length), fill_value=np.NaN)

for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        file = '/Users/jmh2g09/Documents/PhD/Data/Gridded/' + year + month + '_grid.nc'
        nc = Dataset(file, 'r')
        latitude = nc.variables['lat'][:]
        longitude = nc.variables['lon'][:]
        dot = nc.variables['filtered_dynamic_ocean_topography_seasonal_offset'][:]
        nc.close()
        
        nc_grid = Dataset('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/grid.nc', 'w')
        nc_grid.createDimension('lat', np.size(latitude))
        nc_grid.createDimension('lon', np.size(longitude))
        latitudes = nc.createVariable('lat', float, ('lat',))
        longitudes = nc.createVariable('lon', float, ('lon',))
        gridded_dot = nc.createVariable('z', float, ('lat','lon'))
        latitudes[:] = latitude
        longitudes[:] = longitude
        gridded_dot[:] = dot
        nc_grid.close()
        
        f = open('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/coastal_coordinates.dat', 'r')
        ix = 0
        coast_distance = 0
        for line in f:
            line.strip()
            columns = line.split()

            B_lon = columns[0]
            B_lat = columns[1]
            E_lon = columns[2]
            E_lat = columns[3]
            
            # Convert lat/lon to map projection coordinates
            B_x, B_y = m(B_lon, B_lat)
            E_x, E_y = m(E_lon, E_lat)
            
            # Make the line end coordinates into linestring object
            perp_line = LineString(((B_x, B_y), (E_x, E_y)))
            
            # Get interpolated
            line_x = []
            line_y = []
            for iinterp in range(11):
                interp = iinterp/10
                line = perp_line.interpolate(interp, normalized=True).xy
                line_x.append(tuple(line)[0][0])
                line_y.append(tuple(line)[1][0])
            
            # Convert projection x,y coords to lon/lat
            line_lon, line_lat = m(line_x, line_y, inverse=True)
            
            # Convert -180:180 to 0:360
            line_lon = np.array(line_lon)
            line_lat = np.array(line_lat)
            line_lon[np.array(line_lon) < 0] += 360
            
            # Save lon/lat pairs to a file for GMT interpolating 
            input_file = open('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/INPUT.dat', 'w')
            for it in range(len(line_lon)):
                print(line_lon[it], line_lat[it], file=input_file)
            input_file.close()

            # GMT interpolation of dot at the points 
            os.system('gmt grdtrack /Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/INPUT.dat -f0x,1y -G/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/grid.nc > /Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/OUTPUT.dat')
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/INPUT.dat')
            
            # Open the interpolated dot point file
            out_lon = []
            out_lat = []
            out_dot = []
            output_file = open('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/OUTPUT.dat', 'r')
            for line in output_file:
                line.strip()
                columns = line.split()
                if float(columns[0]) < 0:
                    out_lon.append(float(columns[0]) + 360)
                else:
                    out_lon.append(float(columns[0]))
                out_lat.append(float(columns[1]))
                out_dot.append(float(columns[2]))
            output_file.close()
            
            # Save the interpolated dot in an array
            dot_array[int(year)-2011, int(month)-1, ix, :] = np.array(out_dot)
            # Calculate distance along the coastline
            if ix == 0:
                dist = 0
            else:
                dist = vincenty((B_lat, B_lon), (prev_lat, prev_lon)).km
            prev_lat = B_lat
            prev_lon = B_lon
            
            coast_distance += dist
            coast_distance_array[int(year)-2011, int(month)-1, ix] = coast_distance
            
            # Calculate distance off coastline
            for point in range(11):
                if point == 0:
                    dist_2[int(year)-2011, int(month)-1, ix, point] = 0
                else:
                    dist_2[int(year)-2011, int(month)-1, ix, point] = vincenty((out_lat[point], out_lon[point]), (out_lat[0], out_lon[0])).km
            
            ix += 1
        f.close()
        
weddell_start = 0
weddell_end = 36
indian_start = 36
indian_end = 125
ross_start = 125
ross_end = 155
ambel_start = 155
ambel_end = -1

dist_2_mean = np.nanmean(np.nanmean(np.nanmean(dist_2, axis=2), axis=1), axis=0)
coast_distance_mean = np.nanmean(np.nanmean(coast_distance_array, axis=1), axis=0)

dot_zonal_mean = np.nanmean(np.nanmean(np.nanmean(dot_array, axis=2), axis=1), axis=0)

pl.figure()
pl.clf()
pl.plot(dist_2_mean, dot_zonal_mean)
pl.ylabel('DOT (m)')
pl.xlabel('Distance from coast (km)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_zonal_mean.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

dot_perp_mean = np.nanmean(np.nanmean(dot_array, axis=1), axis=0)

pl.figure(figsize=(20, 2))
pl.clf()
pl.pcolor(coast_distance_mean, dist_2_mean, np.ma.masked_invalid(np.transpose(dot_perp_mean)), cmap='RdBu_r')
c = pl.colorbar()
pl.plot([coast_distance_mean[weddell_end], coast_distance_mean[weddell_end]], [dist_2_mean[0], dist_2_mean[-1]], color='k', linestyle='--', linewidth=2)
pl.plot([coast_distance_mean[indian_end], coast_distance_mean[indian_end]], [dist_2_mean[0], dist_2_mean[-1]], color='k', linestyle='--', linewidth=2)
pl.plot([coast_distance_mean[ross_end], coast_distance_mean[ross_end]], [dist_2_mean[0], dist_2_mean[-1]], color='k', linestyle='--', linewidth=2)
pl.plot([coast_distance_mean[ross_end], coast_distance_mean[ross_end]], [dist_2_mean[0], dist_2_mean[-1]], color='k', linestyle='--', linewidth=2)
pl.clim(-1.8, -2.5)
c.set_label('DOT (m)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_mean.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

dy = np.gradient(dist_2, axis=3)

ddot_dy = np.gradient(dot_array, axis=3) / dy

ddot_dy_monthly_mean = np.nanmean(ddot_dy * 1000, axis=0)

pl.figure()
pl.clf()
pl.plot(range(1, 13), np.nanmean(np.nanmean(ddot_dy_monthly_mean[:, :, :4], axis=2), axis=1))
pl.xlim(0, 13)
pl.ylabel('dDOT/do (mm/km)')
pl.xlabel('Month')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_seasonal_gradient.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

####### Split the data into sectors #######

dot_array_monthly_mean = np.nanmean(dot_array, axis=0)

dot_array_weddell = dot_array_monthly_mean[:, weddell_start:weddell_end, :]
ddot_dy_array_weddell = ddot_dy_monthly_mean[:, weddell_start:weddell_end, :]

pl.figure()
pl.clf()
pl.plot(range(1, 13), np.nanmean(np.nanmean(ddot_dy_array_weddell[:, :, :4], axis=2), axis=1))
pl.xlim(0, 13)
pl.ylabel('dDOT/do (mm/km)')
pl.xlabel('Month')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_seasonal_gradient_wedd.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
pl.plot(dist_2_mean, np.nanmean(np.nanmean(dot_array_weddell, axis=1), axis=0))
pl.ylabel('DOT (m)')
pl.xlabel('Distance from coast (km)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_zonal_mean_wedd.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

dot_array_indian = dot_array_monthly_mean[:, indian_start:indian_end, :]
ddot_dy_array_indian = ddot_dy_monthly_mean[:, indian_start:indian_end, :]

pl.figure()
pl.clf()
pl.plot(range(1, 13), np.nanmean(np.nanmean(ddot_dy_array_indian[:, :, :4], axis=2), axis=1))
pl.xlim(0, 13)
pl.ylabel('dDOT/do (mm/km)')
pl.xlabel('Month')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_seasonal_gradient_ind.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
pl.plot(dist_2_mean, np.nanmean(np.nanmean(dot_array_indian, axis=1), axis=0))
pl.ylabel('DOT (m)')
pl.xlabel('Distance from coast (km)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_zonal_mean_ind.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

dot_array_ross = dot_array_monthly_mean[:, ross_start:ross_end, :]
ddot_dy_array_ross = ddot_dy_monthly_mean[:, ross_start:ross_end, :]

pl.figure()
pl.clf()
pl.plot(range(1, 13), np.nanmean(np.nanmean(ddot_dy_array_ross[:, :, :4], axis=2), axis=1))
pl.xlim(0, 13)
pl.ylabel('dDOT/do (mm/km)')
pl.xlabel('Month')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_seasonal_gradient_ross.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
pl.plot(dist_2_mean, np.nanmean(np.nanmean(dot_array_ross, axis=1), axis=0))
pl.ylabel('DOT (m)')
pl.xlabel('Distance from coast (km)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_zonal_mean_ross.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

dot_array_ambel = dot_array_monthly_mean[:, ambel_start:ambel_end, :]
ddot_dy_array_ambel = ddot_dy_monthly_mean[:, ambel_start:ambel_end, :]

pl.figure()
pl.clf()
pl.plot(range(1, 13), np.nanmean(np.nanmean(ddot_dy_array_ambel[:, :, :4], axis=2), axis=1))
pl.xlim(0, 13)
pl.ylabel('dDOT/do (mm/km)')
pl.xlabel('Month')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_seasonal_gradient_ambel.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

pl.figure()
pl.clf()
pl.plot(dist_2_mean, np.nanmean(np.nanmean(dot_array_ambel, axis=1), axis=0))
pl.ylabel('DOT (m)')
pl.xlabel('Distance from coast (km)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_zonal_mean_ambel.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

dot_array_faster = ddot_dy_monthly_mean[4:7, :, :]
dot_array_slower_1 = ddot_dy_monthly_mean[:2, :, :]
dot_array_slower_2 = ddot_dy_monthly_mean[11:, :, :]
dot_array_slower = np.concatenate((dot_array_slower_1, dot_array_slower_2), axis=0)

dot_array_faster_mean = np.nanmean(dot_array_faster, axis=0)
dot_array_slower_mean = np.nanmean(dot_array_slower, axis=0)

dot_array_diff_mean = dot_array_faster_mean - dot_array_slower_mean

pl.figure(figsize=(20, 2))
pl.clf()
pl.pcolor(coast_distance_mean, dist_2_mean, np.ma.masked_invalid(np.transpose(dot_array_diff_mean)), cmap='RdBu_r')
c = pl.colorbar()
pl.clim(0.4, -0.4)
pl.plot([coast_distance_mean[weddell_end], coast_distance_mean[weddell_end]], [dist_2_mean[0], dist_2_mean[-1]], color='k', linestyle='--', linewidth=2)
pl.plot([coast_distance_mean[indian_end], coast_distance_mean[indian_end]], [dist_2_mean[0], dist_2_mean[-1]], color='k', linestyle='--', linewidth=2)
pl.plot([coast_distance_mean[ross_end], coast_distance_mean[ross_end]], [dist_2_mean[0], dist_2_mean[-1]], color='k', linestyle='--', linewidth=2)
pl.plot([coast_distance_mean[ross_end], coast_distance_mean[ross_end]], [dist_2_mean[0], dist_2_mean[-1]], color='k', linestyle='--', linewidth=2)
c.set_label('mm/km')
pl.ylabel('Distance off-coast (km)')
pl.xlabel('Distance along-coast (km)')
pl.title('$\partial$DOT/$\partial$o anomaly (AMJ - NDJF)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/CoastalCurrent/Figures/perp_diff.png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()
