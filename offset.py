import os
import functions as funct
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.basemap import Basemap

yr = input('What year? (xxxx): ')

monthly_offset = []

for mnth in range(1,13):
    offset = []
    if 0 < mnth < 10:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + '0' + str(mnth) + '_elev')
    else:
        os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + str(mnth) + '_elev')

    for file in os.listdir():

        ice_edge = [1] * 10 + [2] * 10
        ice_edge2 = [2] * 10 + [1] * 10

        ssh = []
        surface_type = []

        f = open(file, 'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            # If data point is from open ocean (1) or from a lead (2)
            if columns[0] == '1' or columns[0] == '2':
                # If data point is listed as 'valid'
                if columns[1] == '1':
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
                    offset.append(np.mean(ssh[step - 10:step]) - np.mean(ssh[step:step + 10]))
            if surface_type[step - 10:step] == [2] * 10:
                # Ice to Ocean step
                if abs(np.mean(ssh[step:step + 10]) - np.mean(ssh[step - 10:step])) < 5.:
                    offset.append(np.mean(ssh[step:step + 10]) - np.mean(ssh[step - 10:step]))
        
    monthly_offset.append(np.mean(offset))


print(monthly_offset)

pl.figure()
pl.plot(range(1,13), monthly_offset)
pl.title(yr + ' monthly offset between lead and ocean ssh (m)')
pl.ylabel('ssh offset (m)')
pl.xlabel('Month')
pl.xlim([1, 12])
pl.show()