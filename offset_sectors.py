import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

yr = input('What year? (xxxx): ')

monthly_offset_WEDD = []
monthly_offset_IND = []
monthly_offset_ROSS = []
monthly_offset_AMBEL = []
monthly_offset = []

for mnth in range(1,13):
    offset_WEDD = []
    offset_IND = []
    offset_ROSS = []
    offset_AMBEL = []
    
    if 0 < mnth < 10:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + '0' + str(mnth) + '_elev')
    else:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + str(mnth) + '_elev')

    for file in os.listdir():

        ice_edge = [1] * 10 + [2] * 10
        ice_edge2 = [2] * 10 + [1] * 10

        ssh = []
        lon = []
        surface_type = []

        f = open(file, 'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            # If data point is from open ocean (1) or from a lead (2)
            if columns[0] == '1' or columns[0] == '2':
                # If data point is listed as 'valid'
                if columns[1] == '1':
                    lon.append(float(columns[6]))
                    ssh.append(float(columns[7]))
                    surface_type.append(float(columns[0]))

        iedge = []
        for it in range(len(surface_type)):
            if surface_type[it:it + 20] == ice_edge:
                iedge.append(it + 10)
            if surface_type[it:it + 20] == ice_edge2:
                iedge.append(it + 10)

        for step in iedge:
            if surface_type[step - 10:step] == [1] * 10:
                # Ocean to Ice step
                if abs(np.mean(ssh[step - 10:step]) - np.mean(ssh[step:step + 10])) < 5.:
                    # Choose what the sector the offset point lies within
                    if -60. <= np.mean(lon[step - 10:step + 10]) <= 0.:
                        #print('Weddell')
                        offset_WEDD.append(np.mean(ssh[step - 10:step]) - np.mean(ssh[step:step + 10]))
                    if 0 <= np.mean(lon[step - 10:step + 10]) <= 160.:
                        #print('Indian')
                        offset_IND.append(np.mean(ssh[step - 10:step]) - np.mean(ssh[step:step + 10]))
                    if 160. <= np.mean(lon[step - 10:step + 10]) <= 180.:
                        #print('Ross')
                        offset_ROSS.append(np.mean(ssh[step - 10:step]) - np.mean(ssh[step:step + 10]))
                    if -180. <= np.mean(lon[step - 10:step + 10]) <= -130.:
                        #print('Ross')
                        offset_ROSS.append(np.mean(ssh[step - 10:step]) - np.mean(ssh[step:step + 10]))
                    if -130. <= np.mean(lon[step - 10:step + 10]) <= -60.:
                        #print('Amundsen-Bellingshausen')
                        offset_AMBEL.append(np.mean(ssh[step - 10:step]) - np.mean(ssh[step:step + 10]))
                        
                    
            if surface_type[step - 10:step] == [2] * 10:
                # Ice to Ocean step
                if abs(np.mean(ssh[step:step + 10]) - np.mean(ssh[step - 10:step])) < 5.:
                    # Choose what the sector the offset point lies within
                    if -60. <= np.mean(lon[step - 10:step + 10]) <= 0.:
                        #print('Weddell')
                        offset_WEDD.append(np.mean(ssh[step:step + 10]) - np.mean(ssh[step - 10:step]))
                    if 0 <= np.mean(lon[step - 10:step + 10]) <= 160.:
                        #print('Indian')
                        offset_IND.append(np.mean(ssh[step:step + 10]) - np.mean(ssh[step - 10:step]))
                    if 160. <= np.mean(lon[step - 10:step + 10]) <= 180.:
                        #print('Ross')
                        offset_ROSS.append(np.mean(ssh[step:step + 10]) - np.mean(ssh[step - 10:step]))
                    if -180. <= np.mean(lon[step - 10:step + 10]) <= -130.:
                        #print('Ross')
                        offset_ROSS.append(np.mean(ssh[step:step + 10]) - np.mean(ssh[step - 10:step]))
                    if -130. <= np.mean(lon[step - 10:step + 10]) <= -60.:
                        #print('Amundsen-Bellingshausen')
                        offset_AMBEL.append(np.mean(ssh[step:step + 10]) - np.mean(ssh[step - 10:step]))
        
    monthly_offset_WEDD.append(np.mean(offset_WEDD))
    monthly_offset_IND.append(np.mean(offset_IND))
    monthly_offset_ROSS.append(np.mean(offset_ROSS))  
    monthly_offset_AMBEL.append(np.mean(offset_AMBEL))
    monthly_offset.append((np.mean(offset_WEDD) + np.mean(offset_IND)
        + np.mean(offset_ROSS) + np.mean(offset_AMBEL))/4)

pl.figure()
pl.hist([offset_WEDD, offset_IND, offset_ROSS, offset_AMBEL], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
pl.legend()
pl.title(yr + ' offset histogram')
pl.xlabel('Offset Bin (m)')
pl.ylabel('Frequency')
pl.show()

pl.figure()
pl.plot(range(1,13), monthly_offset, label='Circumpolar')
pl.plot(range(1,13), monthly_offset_WEDD, label='Weddell')
pl.plot(range(1,13), monthly_offset_IND, label='Indian')
pl.plot(range(1,13), monthly_offset_ROSS, label='Ross')
pl.plot(range(1,13), monthly_offset_AMBEL, label='Amundsen-Bellingshausen')
pl.legend(loc='lower right')
pl.title(yr + ' monthly offset between lead and ocean ssh (m)')
pl.ylabel('ssh offset (m)')
pl.xlabel('Month')
pl.xlim([1, 12])
pl.show()