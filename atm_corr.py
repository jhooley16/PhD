import os
import numpy as np
import functions as funct

ionosphere  = {'min': [], 'max': [], 'mean': [], 'std': []}
dry_troposphere  = {'min': [], 'max': [], 'mean': [], 'std': []}
wet_troposphere  = {'min': [], 'max': [], 'mean': [], 'std': []}
dynamic_atmosphere = {'min': [], 'max': [], 'mean': [], 'std': []}
inverse_barometer = {'min': [], 'max': [], 'mean': [], 'std': []}
ocean_tide = {'min': [], 'max': [], 'mean': [], 'std': []}
long_period_tide = {'min': [], 'max': [], 'mean': [], 'std': []}
loading_tide = {'min': [], 'max': [], 'mean': [], 'std': []}
earth_tide = {'min': [], 'max': [], 'mean': [], 'std': []}
pole_tide = {'min': [], 'max': [], 'mean': [], 'std': []}
n = []
n_dyn_atm = []
n_inv_baro = []

for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        directory = '/Volumes/My Passport/Data/atm_cor_files/' + year + month + '_MERGE/'
        print(year, month)
        
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
        
        for file in os.listdir(directory):
            if file[-5:] == '.elev':
                print(file)
                lat = []
                lon = []
                pre_dyn_atm_inv_baro_data = []

                f = open(directory + file, 'r')
                for line in f:
                    line = line.strip()
                    columns = line.split()
                    if len(columns) == 19:
                        # Is the data valid?
                        if columns[1] == '1':
                            # If surface is 'ocean' or 'lead'
                            if (columns[0] == '1') or (columns[0] == '2'):
                                lat.append(float(columns[5]))
                                lon.append(float(columns[6]))
                            
                                ionosphere_data.append(float(columns[10]))
                                dry_troposphere_data.append(float(columns[11]))
                                wet_troposphere_data.append(float(columns[12]))
                            
                                pre_dyn_atm_inv_baro_data.append(float(columns[13]))
                            
                                ocean_tide_data.append(float(columns[14]))
                                long_period_tide_data.append(float(columns[15]))
                                loading_tide_data.append(float(columns[16]))
                                earth_tide_data.append(float(columns[17]))
                                pole_tide_data.append(float(columns[18]))
                f.close()

                # Split the data into inverse barometer and dynamic atmosphere
                # depends on what satellite mode the point was measured with
                pre_dyn_atm_inv_baro_data = np.array(pre_dyn_atm_inv_baro_data)
                mode = np.array(funct.mode_points(lat, lon, month))
            
                dyn_atm_mode = np.where(mode == 0)
                inv_baro_mode = np.where(mode != 0)
            
                dynamic_atmosphere_data = dynamic_atmosphere_data + list(pre_dyn_atm_inv_baro_data[dyn_atm_mode])
                inverse_barometer_data = inverse_barometer_data + list(pre_dyn_atm_inv_baro_data[inv_baro_mode])
                
        # IONOSPHERE
        ionosphere['min'].append(np.round(np.nanmin(ionosphere_data), 3))
        ionosphere['max'].append(np.round(np.nanmax(ionosphere_data), 3))
        ionosphere['mean'].append(np.nanmean(ionosphere_data))
        ionosphere['std'].append(np.std(ionosphere_data))
        # DRY TROPOSPHERE      
        dry_troposphere['min'].append(np.round(np.nanmin(dry_troposphere_data), 3))
        dry_troposphere['max'].append(np.round(np.nanmax(dry_troposphere_data), 3))
        dry_troposphere['mean'].append(np.nanmean(dry_troposphere_data))
        dry_troposphere['std'].append(np.std(dry_troposphere_data))
        # WET TROPOSPHERE
        wet_troposphere['min'].append(np.round(np.nanmin(wet_troposphere_data), 3))
        wet_troposphere['max'].append(np.round(np.nanmax(wet_troposphere_data), 3))
        wet_troposphere['mean'].append(np.nanmean(wet_troposphere_data))
        wet_troposphere['std'].append(np.std(wet_troposphere_data))
        # DYNAMIC ATMOSPHERE
        dynamic_atmosphere['min'].append(np.round(np.nanmin(dynamic_atmosphere_data), 3))
        dynamic_atmosphere['max'].append(np.round(np.nanmax(dynamic_atmosphere_data), 3))
        dynamic_atmosphere['mean'].append(np.nanmean(dynamic_atmosphere_data))
        dynamic_atmosphere['std'].append(np.std(dynamic_atmosphere_data))
        # INVERSE BAROMETER
        inverse_barometer['min'].append(np.round(np.nanmin(inverse_barometer_data), 3))
        inverse_barometer['max'].append(np.round(np.nanmax(inverse_barometer_data), 3))
        inverse_barometer['mean'].append(np.nanmean(inverse_barometer_data))
        inverse_barometer['std'].append(np.std(inverse_barometer_data))
        # OCEAN TIDE
        ocean_tide['min'].append(np.round(np.nanmin(ocean_tide_data), 3))
        ocean_tide['max'].append(np.round(np.nanmax(ocean_tide_data), 3))
        ocean_tide['mean'].append(np.nanmean(ocean_tide_data))
        ocean_tide['std'].append(np.std(ocean_tide_data))
        # LONG TIDE
        long_period_tide['min'].append(np.round(np.nanmin(long_period_tide_data), 3))
        long_period_tide['max'].append(np.round(np.nanmax(long_period_tide_data), 3))
        long_period_tide['mean'].append(np.nanmean(long_period_tide_data))
        long_period_tide['std'].append(np.std(long_period_tide_data))
        # LOADING TIDE
        loading_tide['min'].append(np.round(np.nanmin(loading_tide_data), 3))
        loading_tide['max'].append(np.round(np.nanmax(loading_tide_data), 3))
        loading_tide['mean'].append(np.nanmean(loading_tide_data))
        loading_tide['std'].append(np.std(loading_tide_data))
        # EARTH TIDE
        earth_tide['min'].append(np.round(np.nanmin(earth_tide_data), 3))
        earth_tide['max'].append(np.round(np.nanmax(earth_tide_data), 3))
        earth_tide['mean'].append(np.nanmean(earth_tide_data))
        earth_tide['std'].append(np.std(earth_tide_data))
        # POLE TIDE
        pole_tide['min'].append(np.round(np.nanmin(pole_tide_data), 3))
        pole_tide['max'].append(np.round(np.nanmax(pole_tide_data), 3))
        pole_tide['mean'].append(np.nanmean(pole_tide_data))
        pole_tide['std'].append(np.std(pole_tide_data))
        # NUMBER
        n.append(len(pole_tide_data))
        n_dyn_atm.append(len(dynamic_atmosphere_data))
        n_inv_baro.append(len(inverse_barometer_data))

N = np.sum(n) # Total number of data points
ionosphere_mean = np.nansum(np.array(ionosphere['mean']) * n) / N
ionosphere_std = np.sqrt(np.nansum(n *((np.array(ionosphere['std'])**2) + (np.array(ionosphere['mean']) - ionosphere_mean)**2)) / N)

dry_troposphere_mean = np.nansum(np.array(dry_troposphere['mean']) * n) / N
dry_troposphere_std = np.sqrt(np.nansum(n *((np.array(dry_troposphere['std'])**2) + (np.array(dry_troposphere['mean']) - dry_troposphere_mean)**2)) / N)

wet_troposphere_mean = np.nansum(np.array(wet_troposphere['mean']) * n) / N
wet_troposphere_std = np.sqrt(np.nansum(n *((np.array(wet_troposphere['std'])**2) + (np.array(wet_troposphere['mean']) - wet_troposphere_mean)**2)) / N)

N_dyn_atm = np.sum(n_dyn_atm) # Number of data points (dynamic atmosphere)
dynamic_atmosphere_mean = np.nansum(np.array(dynamic_atmosphere['mean']) * n_dyn_atm) / N_dyn_atm
dynamic_atmosphere_std = np.sqrt(np.nansum(n_dyn_atm *((np.array(dynamic_atmosphere['std'])**2) + (np.array(dynamic_atmosphere['mean']) - dynamic_atmosphere_mean)**2)) / N_dyn_atm)

N_inv_baro = np.sum(n_inv_baro) # Number of data points (inverse barometer)
inverse_barometer_mean = np.nansum(np.array(inverse_barometer['mean']) * n_inv_baro) / N_inv_baro
inverse_barometer_std = np.sqrt(np.nansum(n_inv_baro *((np.array(inverse_barometer['std'])**2) + (np.array(inverse_barometer['mean']) - inverse_barometer_mean)**2)) / N_inv_baro)

ocean_tide_mean = np.nansum(np.array(ocean_tide['mean']) * n) / N
ocean_tide_std = np.sqrt(np.nansum(n *((np.array(ocean_tide['std'])**2) + (np.array(ocean_tide['mean']) - ocean_tide_mean)**2)) / N)

long_period_tide_mean = np.nansum(np.array(long_period_tide['mean']) * n) / N
long_period_tide_std = np.sqrt(np.nansum(n *((np.array(long_period_tide['std'])**2) + (np.array(long_period_tide['mean']) - long_period_tide_mean)**2)) / N)

loading_tide_mean = np.nansum(np.array(loading_tide['mean']) * n) / N
loading_tide_std = np.sqrt(np.nansum(n *((np.array(loading_tide['std'])**2) + (np.array(loading_tide['mean']) - loading_tide_mean)**2)) / N)

earth_tide_mean = np.nansum(np.array(earth_tide['mean']) * n) / N
earth_tide_std = np.sqrt(np.nansum(n *((np.array(earth_tide['std'])**2) + (np.array(earth_tide['mean']) - earth_tide_mean)**2)) / N)

pole_tide_mean = np.nansum(np.array(pole_tide['mean']) * n) / N
pole_tide_std = np.sqrt(np.nansum(n *((np.array(pole_tide['std'])**2) + (np.array(pole_tide['mean']) - pole_tide_mean)**2)) / N)

f = open('/Users/jmh2g09/Documents/PhD/Data/AtmCorrections/atm_corr.dat', 'w')
print('Correction', 'min', 'max', 'mean', 'std', file=f, sep='\t')
print('Ionosphere', np.nanmin(ionosphere['min']), np.nanmax(ionosphere['max']), np.round(ionosphere_mean, 3), np.round(ionosphere_std, 3), file=f, sep='\t')
print('Dry Trop', np.nanmin(dry_troposphere['min']), np.nanmax(dry_troposphere['max']), np.round(dry_troposphere_mean, 3), np.round(dry_troposphere_std, 3), file=f, sep='\t')
print('Wet Trop', np.nanmin(wet_troposphere['min']), np.nanmax(wet_troposphere['max']), np.round(wet_troposphere_mean, 3), np.round(wet_troposphere_std, 3), file=f, sep='\t')
print('Dynamic Atm', np.nanmin(dynamic_atmosphere['min']), np.nanmax(dynamic_atmosphere['max']), np.round(dynamic_atmosphere_mean, 3), np.round(dynamic_atmosphere_std, 3), file=f, sep='\t')
print('Inverse Baro', np.nanmin(inverse_barometer['min']), np.nanmax(inverse_barometer['max']), np.round(inverse_barometer_mean, 3), np.round(inverse_barometer_std, 3), file=f, sep='\t')
print('Ocean Tide', np.nanmin(ocean_tide['min']), np.nanmax(ocean_tide['max']), np.round(ocean_tide_mean, 3), np.round(ocean_tide_std, 3), file=f, sep='\t')
print('Long Tide', np.nanmin(long_period_tide['min']), np.nanmax(long_period_tide['max']), np.round(long_period_tide_mean, 3), np.round(long_period_tide_std, 3), file=f, sep='\t')
print('Loading Tide', np.nanmin(loading_tide['min']), np.nanmax(loading_tide['max']), np.round(loading_tide_mean, 3), np.round(loading_tide_std, 3), file=f, sep='\t')
print('Earth Tide', np.nanmin(earth_tide['min']), np.nanmax(earth_tide['max']), np.round(earth_tide_mean, 3), np.round(earth_tide_std, 3), file=f, sep='\t')
print('Pole Tide', np.nanmin(pole_tide['min']), np.nanmax(pole_tide['max']), np.round(pole_tide_mean, 3), np.round(pole_tide_std, 3), file=f, sep='\t')
f.close()