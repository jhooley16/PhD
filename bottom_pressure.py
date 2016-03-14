import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl

yr = input('What year? (xxxx): ')

temp_day = 0
mdt = []
day_mdt = []

for mnth in range(1, 13):

    if 0 < mnth < 10:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + '0' + str(mnth) + '_elev')
    if 10 <= mnth:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + str(mnth) + '_elev')
    directory = os.getcwd()

    month = directory[-7:-5]
    year = directory[-11:-7]
    print(month, year)
    
    max_day = os.listdir()[-1][13:15]
    for iday in range(1, int(max_day) + 1):
        data = funct.day_data(iday, directory)
        
        lat = data['lat']
        lon = data['lon']
        ssh = data['ssh']
        
        ### If there is no data in the area, assign a NaN to this ###
        if np.size(lat) == 0:
            temp_day = temp_day + 1
        else:
            new_dir = '/Users/jmh2g09/Documents/PhD/Data/BPR/daily_data/'
        
            ### Generate the track File for use with GUT ###
            if 0 < iday < 10:
                nc = Dataset(new_dir + year + month + '0' + str(iday) + '_track.nc', 'w', format='NETCDF3_CLASSIC')
            if 10 <= iday:
                nc = Dataset(new_dir + year + month + str(iday) + '_track.nc', 'w', format='NETCDF3_CLASSIC')

            nc.createDimension('station', np.size(lat))

            longitudes = nc.createVariable('lon', float, ('station',))
            latitudes = nc.createVariable('lat', float, ('station',))
            crs = nc.createVariable('crs', 'i', ())
        
            latitudes.long_name = 'latitude'
            latitudes.standard_name = 'latitude'
            latitudes.units = 'degrees_north'
            longitudes.long_name = 'longitude'
            longitudes.standard_name = 'longitude'
            longitudes.units = 'degrees_east'
            crs.semi_major_axis = 6378137.
            crs.inverse_flattening = 298.257222101004
            crs.earth_gravity_constant = 398600500000000.
            crs.earth_rotation_rate = 7.292115e-05
        
            latitudes[:] = lat
            longitudes[:] = lon
            nc.close()

            if 0 < iday < 10:
                os.system('gut geoidheight_tf -InFile GOCO05s.gfc -T tide-free -InTrack '
                + new_dir + year + month + '0' + str(iday) + '_track.nc -OutFile ' +new_dir 
                + year + month + '0' + str(iday) + '_geoid_track.nc')
            if 10 <= iday:
                os.system('gut geoidheight_tf -InFile GOCO05s.gfc -T tide-free -InTrack '
                + new_dir + year + month + str(iday) + '_track.nc -OutFile ' + new_dir +
                year + month + str(iday) + '_geoid_track.nc')
        
            ### Open the geoid_track file and compute the mdt ###
        
            if 0 < iday < 10:
                nc = Dataset(new_dir + year + month + '0' + str(iday) + '_geoid_track.nc', 'r')
            if 10 <= iday:
                nc = Dataset(new_dir + year + month + str(iday) + '_geoid_track.nc', 'r')
        
            geoid_height = nc.variables['geoid_height'][:]
        
            nc.close()
        
            mdt.append(np.nanmean(ssh - geoid_height))
            temp_day = temp_day + 1
            day_mdt.append(temp_day)

mdta = mdt - np.nanmean(mdt)

# Bottom pressure recorder located in the Weddell Sea
# -60.8249 N
# -54.7221 E (305.2779 W)
# depth = 1920 m 
file = '/Users/jmh2g09/Documents/PhD/Data/BPR/DPS_DEEP_1113_DQ105443_drp.txt'

day = []
bpa_millibars = []

f = open(file, 'r')

for line in f:
    line = line.strip()
    columns = line.split()
    if np.size(columns) > 4:
        # only choose data for the year 2012
        if columns[2] == '2012':
            # only choose 'valid' data
            if columns[1] == '0':
                day.append(float(columns[3]))
                bpa_millibars.append(float(columns[5]))

# bottom pressure is measured in millibars, need to convert to Pascals
# 1 millibar = 100 Pa
bpa_pascals = np.array(bpa_millibars) * 100

# Calculate the approximate change in sea level to produce the pressure change
dh = bpa_pascals / (1025 * 9.81)

match_dh = []
for iday in day_mdt:
    if iday == day:
        I = day.index(iday)
        match_dh.append(dh[I])

print('the size of dh: ', np.size(match_dh), 'the size of mtda', np.size(mdta))
pl.figure()
pl.scatter(match_dh, mdta)
pl.show()

# create the general figure
fig1 = pl.figure()
# and the first axes using subplot populated with data 
ax1 = fig1.add_subplot(111)
line1 = ax1.plot(day, funct.smooth(dh, window_len=30), color='b')
pl.ylabel('Bottom Pressure Anomaly (m)')
 
# now, the second axes that shares the x-axis with the ax1
ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
line2 = ax2.plot(day_mdt, mdta, color='r')
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
pl.ylabel('Sea Surface Height Anomaly (m)')
 
# for the legend, remember that we used two different axes so, we need 
# to build the legend manually
pl.legend((line1, line2), ('Bottom Pressure Anomaly (m)', 'Sea Surface Height Anomaly (m)'))
pl.show()