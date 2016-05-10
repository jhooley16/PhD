import os
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

yr = '2015'#input('What year? (xxxx) ')

correlation = []
corr_length = []
distance = []

for mnth in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:

    if os.path.isdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + mnth + '_elev'):
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + mnth + '_elev')
        tracks = os.listdir()
        for file in tracks:
            #print(file)

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

            # Plot the satellite tracks
            #pl.figure()
            #pl.clf()
            #m.drawmapboundary()
            #m.drawcoastlines(zorder=10)
            #m.fillcontinents(zorder=10)
            #m.drawparallels(np.arange(-80., 81., 20.), labels=[1, 0, 0, 0])
            #m.drawmeridians(np.arange(-180., 181., 20.), labels=[0, 0, 0, 1])

            # Calculate the autocorrelation of the ASCENDING track
            if len(ascending) > 0:
                ssha_asc = ssha[ascending[1]:ascending[-1]]
                lat_asc = lat[ascending[1]:ascending[-1]]
                lon_asc = lon[ascending[1]:ascending[-1]] 
                           
                if abs(lat_asc[0] - lat_asc[-1]) > 15.:
                    if np.max(abs(np.gradient(lat_asc))) < .1:
                        stereo_x, stereo_y = m(lon_asc, lat_asc)
                        #m.scatter(stereo_x, stereo_y, color='b')
                
                        d_asc = ((np.array(stereo_x)**2) + (np.array(stereo_y)**2))**0.5
                        distance.append(np.mean(abs(np.gradient(d_asc))))

                        corr_asc = np.correlate(ssha_asc, ssha_asc, mode='same')
                        N_asc = len(corr_asc)
                        half_asc = corr_asc[N_asc//2:]
                        lengths_asc = range(N_asc, N_asc//2, -1)
                        half_asc /= lengths_asc
                        half_asc /= half_asc[0]

                        correlation.append(half_asc)
                        corr_length.append(len(half_asc))

            # Calculate the autocorrelation of the DESCENDING track
            if len(descending) > 0:
                ssha_desc = ssha[descending[0]:descending[-1]]
                lat_desc = lat[descending[0]:descending[-1]]
                lon_desc = lon[descending[0]:descending[-1]]
                
                if abs(lat_desc[0] - lat_desc[-1]) > 15.:
                    if np.max(abs(np.gradient(lat_desc))) < .1:
                        stereo_x, stereo_y = m(lon_desc, lat_desc)
                        #m.scatter(stereo_x, stereo_y, color='r')
                
                        d_desc = ((np.array(stereo_x)**2) + (np.array(stereo_y)**2))**0.5
                        distance.append(np.mean(abs(np.gradient(d_desc))))

                        corr_desc = np.correlate(ssha_desc, ssha_desc, mode='same')
                        N_desc = len(corr_desc)
                        half_desc = corr_desc[N_desc//2:]
                        lengths_desc = range(N_desc, N_desc//2, -1)
                        half_desc /= lengths_desc
                        half_desc /= half_desc[0]
        
                        correlation.append(half_desc)
                        corr_length.append(len(half_desc))

            #pl.show()
            #pl.close()

# Get all the correlations to the same length
A = np.where(corr_length == np.max(corr_length))[0]
sum_correlation = np.zeros(np.shape(correlation[A]))

for icorr in range(0, len(correlation)):
    sum_correlation = sum_correlation + np.resize(correlation[icorr], np.shape(correlation[A]))

mean_correlation = sum_correlation / len(correlation)

print('The mean distance between points is: ' + str(np.mean(distance)/1000) + ' km')
print('The number of orbits used is: ' + str(len(correlation)//2))
pl.figure()
pl.plot(range(len(mean_correlation)) * np.mean(distance)/1000, mean_correlation)
pl.xlim([0, 150])
pl.show()