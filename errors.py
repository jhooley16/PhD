import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
import os
from datetime import date

lat_resolution = '0.5'    #input('What latitude resolution?: ')
lon_resolution = '1.0'      #input ('What longitude resolution?: ')

grid_SSH_RMS_months = np.full((181, 30, 6, 12), fill_value=np.NaN)
grid_SSH_RMS_all = np.full((181, 30, 72), fill_value=np.NaN)
ix = 0
dates = []
for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        dates.append(date(int(year), int(month), 15))
        # Cycle through each raw file
        file = '/Users/jmh2g09/Documents/PhD/Data/Processed/' + year + month + '_track.nc'
        print(year, month)
        nc = Dataset(file, 'r')
        lat = nc.variables['latitude'][:]
        lon = nc.variables['longitude'][:]
        ssh = nc.variables['sea_surface_height'][:]
        time = nc.variables['time'][:]
        mode = nc.variables['mode'][:]
        surface = nc.variables['surface'][:]
        nc.close()

        input_ssh = open('INPUT_ssh.dat', 'w')
        input_t = open('INPUT_t.dat', 'w')
        input_m = open('INPUT_m.dat', 'w')
        input_s = open('INPUT_s.dat', 'w')
        for ilen in range(len(lat)):
            print(lon[ilen], lat[ilen], ssh[ilen], file=input_ssh)
            print(lon[ilen], lat[ilen], time[ilen], file=input_t)
            print(lon[ilen], lat[ilen], mode[ilen], file=input_m)
            print(lon[ilen], lat[ilen], surface[ilen], file=input_s)
        input_ssh.close()
        input_t.close()
        input_m.close()
        input_s.close()

        # Calculate the number of points per node
        os.system('gmt xyz2grd INPUT_ssh.dat -GOUTPUT_n.nc -An -I0.25/0.002 -R-180/180/-79/-50 -fig')

        # Get the maximum and minimum ssh values in each node
        os.system('gmt xyz2grd INPUT_ssh.dat -GOUTPUT_ssh_u.nc -Au -I0.25/0.002 -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_ssh.dat -GOUTPUT_ssh_l.nc -Al -I0.25/0.002 -R-180/180/-79/-50 -fig')
        # Get surface modes
        os.system('gmt xyz2grd INPUT_m.dat -GOUTPUT_m.nc -Au -I0.25/0.002 -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_s.dat -GOUTPUT_s.nc -Au -I0.25/0.002 -R-180/180/-79/-50 -fig')
        # Get minimum and maximum times, to calculate time difference
        os.system('gmt xyz2grd INPUT_t.dat -GOUTPUT_t_u.nc -Au -I0.25/0.002 -R-180/180/-79/-50 -fig')
        os.system('gmt xyz2grd INPUT_t.dat -GOUTPUT_t_l.nc -Al -I0.25/0.002 -R-180/180/-79/-50 -fig')
        os.system('rm INPUT_ssh.dat INPUT_m.dat INPUT_s.dat INPUT_t.dat')
        
        # Load number data
        nc = Dataset('OUTPUT_n.nc', 'r')
        grid_lat = nc.variables['lat'][:]
        grid_lon = nc.variables['lon'][:]
        grid_n = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        # Load minimum and maximum ssh data
        nc = Dataset('OUTPUT_ssh_u.nc', 'r')
        grid_ssh_u = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        nc = Dataset('OUTPUT_ssh_l.nc', 'r')
        grid_ssh_l = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        # Load surface data
        nc = Dataset('OUTPUT_m.nc', 'r')
        grid_m = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        nc = Dataset('OUTPUT_s.nc', 'r')
        grid_s = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        
        # Load time data
        nc = Dataset('OUTPUT_t_u.nc', 'r')
        grid_t_u = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        nc = Dataset('OUTPUT_t_l.nc', 'r')
        grid_t_l = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        os.system('rm OUTPUT_n.nc OUTPUT_ssh_u.nc OUTPUT_ssh_l.nc OUTPUT_m.nc OUTPUT_s.nc OUTPUT_t_u.nc OUTPUT_t_l.nc')
        
        
        ##### Calculate differences #####
        # Take the difference between ssh nodes
        grid_ssh_diff = abs(grid_ssh_u - grid_ssh_l)
        
        # Take the difference between the time nodes
        grid_t_diff = abs(grid_t_u - grid_t_l)
        
        ##### Node Removal #####
        # Remove nodes that do not have 2 data points
        grid_ssh_diff[grid_n != 2] = np.NaN
        
        # Remove nodes that have data points measured less than 10 days apart
        grid_ssh_diff[grid_t_diff >= 10] = np.NaN
        
        # Separate nodes based on surface (open ocean, lead)
        grid_ssh_ocean = grid_ssh_diff[grid_s == 1]
        grid_ssh_lead = grid_ssh_diff[grid_s == 2]
        
        # Separate nodes based on satellite mode (LRM, SAR, SARIn)
        grid_ssh_LRM = grid_ssh_diff[grid_m == 0]
        grid_ssh_SAR = grid_ssh_diff[grid_m == 1]
        grid_ssh_SARIn = grid_ssh_diff[grid_m == 2]
        
        # Convert this grid into an xyz table for re-binning onto regular grid
        grid_lats, grid_lons = np.meshgrid(grid_lat, grid_lon)
        
        input_SSH = open('INPUT_SSH.dat', 'w')
        for ilon in range(np.shape(grid_ssh_diff)[0]):
            for ilat in range(np.shape(grid_ssh_diff)[1]):
                LON = grid_lons[ilon, ilat]
                LAT = grid_lats[ilon, ilat]
                SSH = grid_ssh_diff[ilon, ilat]
                if np.isfinite(SSH):
                    print(LON, LAT, SSH, file=input_SSH)
        input_SSH.close()
        
        # Regrid the data to standard grid
        os.system('gmt xyz2grd INPUT_SSH.dat -GOUTPUT_SSH.nc -Ar -I2.0/1.0 -R-180/180/-79/-50 -fig')
        os.system('rm INPUT_SSH.dat')
        
        # Open regridded data and construct arrays
        nc = Dataset('OUTPUT_SSH.nc', 'r')
        grid_lat_2 = nc.variables['lat'][:]
        grid_lon_2 = nc.variables['lon'][:]
        
        grid_SSH_RMS_months[:, :, int(year)-2011, int(month)-1] = np.array(np.transpose(nc.variables['z'][:]))
        grid_SSH_RMS_all[:, :, ix] = np.array(np.transpose(nc.variables['z'][:]))
        ix += 1
        
        grid_SSH_RMS = np.array(np.transpose(nc.variables['z'][:]))
        nc.close()
        os.system('rm OUTPUT_SSH.nc')

# Take monthly and annual means
grid_SSH_RMS_month = np.nanmean(grid_SSH_RMS_months, 2)
grid_SSH_RMS_years = np.nanmean(grid_SSH_RMS_months, 3)
grid_SSH_RMS_period = np.nanmean(grid_SSH_RMS_years, 2)

for it in range(12):
    
    pl.figure()
    pl.clf()
    m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
    m.drawmapboundary()
    m.drawcoastlines(zorder=10)
    m.fillcontinents(zorder=10)
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
    grid_lats, grid_lons = np.meshgrid(grid_lat_2, grid_lon_2)
    stereo_x, stereo_y = m(grid_lons, grid_lats)
    m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(grid_SSH_RMS_month[:, :, it]))
    m.colorbar()
    pl.clim(0, 0.2)
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Errors/cross_overs_' + str(it) + '.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
    pl.close()


pl.figure()
pl.clf()
m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')
m.drawmapboundary()
m.drawcoastlines(zorder=10)
m.fillcontinents(zorder=10)
m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])
grid_lats, grid_lons = np.meshgrid(grid_lat_2, grid_lon_2)
stereo_x, stereo_y = m(grid_lons, grid_lats)
m.pcolor(stereo_x, stereo_y, np.ma.masked_invalid(grid_SSH_RMS_period))
m.colorbar()
pl.clim(0, 0.2)
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Errors/period_cross_overs.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

# Average and produce regional average time series
period_ts = np.nanmean(np.nanmean(grid_SSH_RMS_all, 1), 0)
period_spread = np.std(np.std(grid_SSH_RMS_all, 1), 0)

fig = pl.figure()
pl.fill_between(dates, period_ts-period_spread, period_ts+period_spread)
pl.plot(dates, period_ts)
fig.autofmt_xdate()
pl.ylabel('Crossover RMS (m)')
pl.xlabel('Date')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Errors/period_ts.png', format='png', transparent=True, dpi=300, bbox_inches='tight')
pl.close()

print('Average RMS crossover difference (2011 - 2016): ' + str(np.nanmean(period_ts)) + ' m')

# Average RMS crossover difference (2011 - 2016): 0.100301367567 m
