import os
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap
from math import exp

correlation = []
corr_length = []
distance = []

for year in ['2011', '2012', '2013', '2014', '2015', '2016']:
    print(year)
    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        print(month)
        directory = '/Volumes/My Passport/Data/elev_files/' + year + month + '_MERGE'
        os.chdir(directory)
        for file in os.listdir():
            ssha = []
            lon = []
            lat = []
            surface_type = []
            time = []

            f = open(file, 'r')
        
            for line in f:
                line = line.strip()
                columns = line.split()
                # If data point is from open ocean (1) or from a lead (2)
                if columns[0] == '1':#or columns[0] == '2':
                    # If data point is listed as 'valid'
                    if columns[1] == '1':
                        if float(columns[7]) - float(columns[8]) < 3.:
                            lat.append(float(columns[5]))
                            lon.append(float(columns[6]))
                            ssha.append(float(columns[7]) - float(columns[8]))
                            surface_type.append(float(columns[0]))

            # Seperate ascending and descending satellite tracks
            descending = np.where(np.gradient(lat) < 0.)[0]
            ascending = np.where(np.gradient(lat) > 0.)[0]
        
            # Calculate the mean difference between satellite track points
            m = Basemap(projection='spstere', boundinglat=-50, lon_0=180, resolution='l')

            # Calculate the autocorrelation of the ASCENDING track
            if len(ascending) > 0:
                ssha_asc = ssha[ascending[1]:ascending[-1]]
                lat_asc = lat[ascending[1]:ascending[-1]]
                lon_asc = lon[ascending[1]:ascending[-1]]
                       
                if np.max(abs(np.gradient(lat_asc))) < .01:
                    if abs(lat_asc[0] - lat_asc[-1]) > 20.:
                        stereo_x, stereo_y = m(lon_asc, lat_asc)
                
                        d_asc = ((np.array(stereo_x)**2) + (np.array(stereo_y)**2))**0.5
                        distance.append(np.mean(abs(np.gradient(d_asc))))

                        corr_asc = np.correlate(ssha_asc, ssha_asc, mode='same')
                        N_asc = len(corr_asc)
                        corr_asc /= corr_asc[N_asc//2]

                        correlation.append(corr_asc)
                        corr_length.append(len(corr_asc))

            # Calculate the autocorrelation of the DESCENDING track
            if len(descending) > 0:
                ssha_desc = ssha[descending[0]:descending[-1]]
                lat_desc = lat[descending[0]:descending[-1]]
                lon_desc = lon[descending[0]:descending[-1]]
            
                if np.max(abs(np.gradient(lat_desc))) < .01:
                    if abs(lat_desc[0] - lat_desc[-1]) > 20.:
                        stereo_x, stereo_y = m(lon_desc, lat_desc)
            
                        d_desc = ((np.array(stereo_x)**2) + (np.array(stereo_y)**2))**0.5
                        distance.append(np.mean(abs(np.gradient(d_desc))))

                        corr_desc = np.correlate(ssha_desc, ssha_desc, mode='same')
                        N_desc = len(corr_desc)
                        corr_desc /= corr_desc[N_desc//2]
    
                        correlation.append(corr_desc)
                        corr_length.append(len(corr_desc))

# Get all the correlations to the same length
A = np.where(corr_length == np.min(corr_length))[0]
sum_correlation = np.zeros(np.shape(correlation[A]))

for icorr in range(0, len(correlation)):
    sum_correlation += correlation[icorr][(len(correlation[icorr])//2) - len(sum_correlation)//2:(len(correlation[icorr])//2) + 1 + len(sum_correlation)//2]

mean_correlation = sum_correlation / len(correlation)

print('The mean distance between points is: ' + str(np.mean(distance)/1000) + ' km')
# 0.252 km
print('The number of orbits used is: ' + str(len(correlation)//2))
# 513
pl.figure()
pl.plot(range(-len(mean_correlation)//2, len(mean_correlation)//2) * np.mean(distance)/1000, mean_correlation)

pl.xlim([-100, 100])
pl.xlabel('lag (km)')
pl.ylabel('Correlation Coefficient')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Offset/Figures/autocorr.png', format='png', transparent=True, doi=300, bbox_inches='tight')
pl.close()