from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
from datetime import date
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import os

GRACE_file = '/Users/jmh2g09/Documents/PhD/Data/GRACE/GRCTellus.JPL.200208_201603.OCN.RL05_1.DSTvDPC1412.nc'

## Load the GRACE data from the raw file. 
nc = Dataset(GRACE_file, 'r')
GRACE_lat = nc.variables['lat'][:]
GRACE_lon = nc.variables['lon'][:]
GRACE_time = nc.variables['time_bounds'][:]
# Convert GRACE lwe thickness to m (from cm) 
GRACE_data_pre = nc.variables['lwe_thickness'][:] / 10 # (time, lat, lon)
nc.close()

# Mask out missing values (land, etc) 
GRACE_data_pre[GRACE_data_pre<-999] = np.NaN

## Load the GRACE data within the time bounds that I am working with
GRACE_year = []
GRACE_month = []
GRACE_day = []
GRACE_data = []
test = False
t_correction = date.toordinal(date(2002, 1, 1))
for i in range(len(GRACE_time)):
    iyear_start = date.timetuple(date.fromordinal(int(GRACE_time[i, 0]) + t_correction))[0]
    imonth_start = date.timetuple(date.fromordinal(int(GRACE_time[i, 0]) + t_correction))[1]
    iday_start = date.timetuple(date.fromordinal(int(GRACE_time[i, 0]) + t_correction))[2]
    iyear_end = date.timetuple(date.fromordinal(int(GRACE_time[i, 1]) + t_correction))[0]
    imonth_end = date.timetuple(date.fromordinal(int(GRACE_time[i, 1]) + t_correction))[1]
    iday_end = date.timetuple(date.fromordinal(int(GRACE_time[i, 1]) + t_correction))[2]
    
    if iyear_start == 2010 and imonth_start == 11:
        test = True
    if iyear_start == 2016 and imonth_start == 3:
        test = False
    
    if test == True:
        GRACE_year.append([iyear_start, iyear_end])
        GRACE_month.append([imonth_start, imonth_end])
        GRACE_day.append([iday_start, iday_end])
        GRACE_data.append(GRACE_data_pre[i, :, :])

GRACE_data = np.array(GRACE_data)

## Load the missing month data
missing_years = []
missing_months = []
missing_file = open('/Users/jmh2g09/Documents/PhD/Data/GRACE/GRACE_missing_months.csv', 'r')
for line in missing_file:
    row = line.strip()
    columns = row.split(',')
    missing_years.append(columns[0])
    missing_months.append(columns[1])

## Find the indices where the missing months are located
grace_ts = []
month_ts = []
trigger = False
month_step = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    # Cycle through the months
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if year == '2010' and month == '11':
            trigger = True
        if year == '2016' and month == '03':
            trigger = False
        
        if trigger == True:
            MISSING = False
            for imiss in range(len(missing_years)):
                if int(year) == int(missing_years[imiss]) and int(month) == int(missing_months[imiss]):
                    MISSING = True

            if MISSING == True:
                pass
            else:
                grace_ts.append(month_step)
            month_ts.append(month_step)
            month_step += 1

## Interpolate the missing months
GRACE = np.full((64, 180, 360), fill_value=np.NaN)
for ilat in range(180):
    for ilon in range(360):
        interp_ts = interpolate.interp1d(grace_ts, GRACE_data[:, ilat, ilon], kind='cubic')
        GRACE[:, ilat, ilon] = interp_ts(month_ts)

## Save the interpolated data in seperate .nc files
## Also resample onto a 0.5 lat, 1.0 lon grid

trigger == False
A = 0
for year in ['2010', '2011', '2012', '2013', '2014', '2015', '2016']:
    # Cycle through the months
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        if year == '2010' and month == '11':
            trigger = True
        if year == '2016' and month == '03':
            trigger = False
        
        if trigger == True:
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/GRACE/INPUT.nc', 'w')
            nc.createDimension('lat', np.size(GRACE_lat))
            nc.createDimension('lon', np.size(GRACE_lon))
            latitude = nc.createVariable('lat', float, ('lat',))
            longitude = nc.createVariable('lon', float, ('lon',))
            GRACE_save = nc.createVariable('GRACE', float, ('lat','lon',))
            latitude[:] = GRACE_lat
            longitude[:] = GRACE_lon
            GRACE_save[:] = GRACE[A, :, :]
            nc.close()
            
            os.system('gmt grdsample /Users/jmh2g09/Documents/PhD/Data/GRACE/INPUT.nc -G/Users/jmh2g09/Documents/PhD/Data/GRACE/OUTPUT.nc -I1.0/0.5 -R0/360/-79/-50')
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/GRACE/INPUT.nc')
            
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/GRACE/OUTPUT.nc', 'r')
            lat = nc.variables['y'][:]
            lon = nc.variables['x'][:]
            z = nc.variables['z'][:]
            nc.close()
            os.system('rm /Users/jmh2g09/Documents/PhD/Data/GRACE/OUTPUT.nc')
            
            nc = Dataset('/Users/jmh2g09/Documents/PhD/Data/GRACE/' + year + month + '_GRACE.nc', 'w')

            nc.createDimension('lat', np.size(lat))
            nc.createDimension('lon', np.size(lon))

            latitude = nc.createVariable('latitude', float, ('lat',))
            longitude = nc.createVariable('longitude', float, ('lon',))
            GRACE_save = nc.createVariable('GRACE', float, ('lat','lon',))

            latitude.long_name = 'latitude'
            latitude.standard_name = 'latitude'
            latitude.units = 'degrees_north'
            longitude.long_name = 'longitude'
            longitude.standard_name = 'longitude'
            longitude.units = 'degrees_east'
            GRACE_save.long_name = 'GRACE measurements'
            GRACE_save.standard_name = 'ocean_gravity_changes_equivalent_water_thickness'
            GRACE_save.units = 'm'

            latitude[:] = lat
            longitude[:] = lon
            GRACE_save[:] = z

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
        
            m.pcolor(stereo_x, stereo_y, np.transpose(np.ma.masked_invalid(z)), cmap='RdBu_r')
            m.colorbar()
            pl.savefig('/Users/jmh2g09/Documents/PhD/Data/GRACE/Figures/' + year + month + '_GRACE.png', transparent=True, dpi=300)
            pl.close()
            
            A += 1


## TODO: Correct for atmospheric changes in pressure, which result in an increase in water depth and therefore mass