import os
import numpy as np

dir = '/Users/jmh2g09/Documents/PhD/Data/AtmCorrections'
os.chdir(dir)

mode = input('What Mode (LRM, SAR, SARIn)?: ')

if mode == 'LRM':
    dir = '/Users/jmh2g09/Documents/PhD/Data/AtmCorrections/201408/201408_LRM'
    os.chdir(dir)
elif mode == 'SAR':
    dir = '/Users/jmh2g09/Documents/PhD/Data/AtmCorrections/201408/201408_SAR'
    os.chdir(dir)
elif mode == 'SARIn':
    dir = '/Users/jmh2g09/Documents/PhD/Data/AtmCorrections/201408/201408_SARIN'
    os.chdir(dir)

files = os.listdir(dir)

ionosphere_data = []
dry_troposphere_data = []
wet_troposphere_data = []
dynamic_atmosphere_data = []
inverse_barometer_data = []
ocean_tide_data = []
long_period_tide_data = []
loading_tide_data = []
earth_tide_data = []
pole_tide_data = []
for file in files:
    if file[-5:] == '.elev':
        f = open(file, 'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            # Is the data valid?
            if columns[1] == '1':
                # If surface is 'ocean' or 'lead'
                if (columns[0] == '1') or (columns[0] == '2'):
                    ionosphere_data.append(float(columns[10]))
                    dry_troposphere_data.append(float(columns[11]))
                    wet_troposphere_data.append(float(columns[12]))
                    if mode == 'LRM':
                        dynamic_atmosphere_data.append(float(columns[13]))
                    else:
                        inverse_barometer_data.append(float(columns[13]))
                    ocean_tide_data.append(float(columns[14]))
                    long_period_tide_data.append(float(columns[15]))
                    loading_tide_data.append(float(columns[16]))
                    earth_tide_data.append(float(columns[17]))
                    pole_tide_data.append(float(columns[18]))

ionosphere = {'min': np.round(np.nanmin(ionosphere_data), 3), 'max': np.round(np.nanmax(ionosphere_data), 3), 'mean': np.round(np.nanmean(ionosphere_data), 3), 'std': np.round(np.std(ionosphere_data), 3)}
dry_troposphere  = {'min': np.round(np.nanmin(dry_troposphere_data), 3), 'max': np.round(np.nanmax(dry_troposphere_data), 3), 'mean': np.round(np.nanmean(dry_troposphere_data), 3), 'std': np.round(np.std(dry_troposphere_data), 3)}
wet_troposphere = {'min': np.round(np.nanmin(wet_troposphere_data), 3), 'max': np.round(np.nanmax(wet_troposphere_data), 3), 'mean': np.round(np.nanmean(wet_troposphere_data), 3), 'std': np.round(np.std(wet_troposphere_data), 3)}
if mode == 'LRM':
    dynamic_atmosphere = {'min': np.round(np.nanmin(dynamic_atmosphere_data), 3), 'max': np.round(np.nanmax(dynamic_atmosphere_data), 3), 'mean': np.round(np.nanmean(dynamic_atmosphere_data), 3), 'std': np.round(np.std(dynamic_atmosphere_data), 3)}
else:
    inverse_barometer = {'min': np.round(np.nanmin(inverse_barometer_data), 3), 'max': np.round(np.nanmax(inverse_barometer_data), 3), 'mean': np.round(np.nanmean(inverse_barometer_data), 3), 'std': np.round(np.std(inverse_barometer_data), 3)}
ocean_tide = {'min': np.round(np.nanmin(ocean_tide_data), 3), 'max': np.round(np.nanmax(ocean_tide_data), 3), 'mean': np.round(np.nanmean(ocean_tide_data), 3), 'std': np.round(np.std(ocean_tide_data), 3)}
long_period_tide = {'min': np.round(np.nanmin(long_period_tide_data), 3), 'max': np.round(np.nanmax(long_period_tide_data), 3), 'mean': np.round(np.nanmean(long_period_tide_data), 3), 'std': np.round(np.std(long_period_tide_data), 3)}
loading_tide = {'min': np.round(np.nanmin(loading_tide_data), 3), 'max': np.round(np.nanmax(loading_tide_data), 3), 'mean': np.round(np.nanmean(loading_tide_data), 3), 'std': np.round(np.std(loading_tide_data), 3)}
earth_tide = {'min': np.round(np.nanmin(earth_tide_data), 3), 'max': np.round(np.nanmax(earth_tide_data), 3), 'mean': np.round(np.nanmean(earth_tide_data), 3), 'std': np.round(np.std(earth_tide_data), 3)}
pole_tide = {'min': np.round(np.nanmin(pole_tide_data), 3), 'max': np.round(np.nanmax(pole_tide_data), 3), 'mean': np.round(np.nanmean(pole_tide_data), 3), 'std': np.round(np.std(pole_tide_data), 3)}

f = open('/Users/jmh2g09/Documents/PhD/Data/AtmCorrections/corr_' + mode + '.dat', 'w')
print('Correction', 'min', 'max', 'mean', 'std', file=f, sep='\t')
print('Ionosphere', ionosphere['min'], ionosphere['max'], ionosphere['mean'], ionosphere['std'], file=f, sep='\t')
print('Dry Trop', dry_troposphere['min'], dry_troposphere['max'], dry_troposphere['mean'], dry_troposphere['std'], file=f, sep='\t')
print('Wet Trop', wet_troposphere['min'], wet_troposphere['max'], wet_troposphere['mean'], wet_troposphere['std'], file=f, sep='\t')
if mode == 'LRM':
    print('Dynamic Atm', dynamic_atmosphere['min'], dynamic_atmosphere['max'], dynamic_atmosphere['mean'], dynamic_atmosphere['std'], file=f, sep='\t')
else:
    print('Inverse Baro', inverse_barometer['min'], inverse_barometer['max'], inverse_barometer['mean'], inverse_barometer['std'], file=f, sep='\t')
print('Ocean Tide', ocean_tide['min'], ocean_tide['max'], ocean_tide['mean'], ocean_tide['std'], file=f, sep='\t')
print('Long Tide', long_period_tide['min'], long_period_tide['max'], long_period_tide['mean'], long_period_tide['std'], file=f, sep='\t')
print('Loading Tide', loading_tide['min'], loading_tide['max'], loading_tide['mean'], loading_tide['std'], file=f, sep='\t')
print('Earth Tide', earth_tide['min'], earth_tide['max'], earth_tide['mean'], earth_tide['std'], file=f, sep='\t')
print('Pole Tide', pole_tide['min'], pole_tide['max'], pole_tide['mean'], pole_tide['std'], file=f, sep='\t')
f.close()
