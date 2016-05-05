import os
import functions as funct
import numpy as np
import matplotlib.pyplot as pl

for year in range(2010, 2017):
    print(year)
    yr = str(year)
    
    monthly_offset_WEDD = []
    monthly_offset_IND = []
    monthly_offset_ROSS = []
    monthly_offset_AMBEL = []
    monthly_offset = []

    hist_offset_WEDD = []
    hist_offset_IND = []
    hist_offset_ROSS = []
    hist_offset_AMBEL = []

    month_number = []

    for mnth in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:

        offset_WEDD = []
        offset_IND = []
        offset_ROSS = []
        offset_AMBEL = []

        if os.path.isdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + mnth + '_elev'):
            os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + mnth + '_elev')
            month_number.append(int(mnth))
        
            for file in os.listdir():

                ice_edge = [1] * 10 + [2] * 10
                ice_edge2 = [2] * 10 + [1] * 10

                ssha = []
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
                            ssha.append(float(columns[7]) - float(columns[8]))
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
                        if abs(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10])) < .5:
                            # Choose what the sector the offset point lies within
                            if -60. <= np.mean(lon[step - 10:step + 10]) <= 0.:
                                #print('Weddell')
                                offset_WEDD.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))
                                hist_offset_WEDD.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))
                            if 0 <= np.mean(lon[step - 10:step + 10]) <= 160.:
                                #print('Indian')
                                offset_IND.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))
                                hist_offset_IND.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))
                            if 160. <= np.mean(lon[step - 10:step + 10]) <= 180.:
                                #print('Ross')
                                offset_ROSS.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))
                                hist_offset_ROSS.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))
                            if -180. <= np.mean(lon[step - 10:step + 10]) <= -130.:
                                #print('Ross')
                                offset_ROSS.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))
                                hist_offset_ROSS.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))
                            if -130. <= np.mean(lon[step - 10:step + 10]) <= -60.:
                                #print('Amundsen-Bellingshausen')
                                offset_AMBEL.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))
                                hist_offset_AMBEL.append(np.mean(ssha[step - 10:step]) - np.mean(ssha[step:step + 10]))


                    if surface_type[step - 10:step] == [2] * 10:
                        # Ice to Ocean step
                        if abs(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step])) < .5:
                            # Choose what the sector the offset point lies within
                            if -60. <= np.mean(lon[step - 10:step + 10]) <= 0.:
                                #print('Weddell')
                                offset_WEDD.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))
                                hist_offset_WEDD.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))
                            if 0 <= np.mean(lon[step - 10:step + 10]) <= 160.:
                                #print('Indian')
                                offset_IND.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))
                                hist_offset_IND.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))
                            if 160. <= np.mean(lon[step - 10:step + 10]) <= 180.:
                                #print('Ross')
                                offset_ROSS.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))
                                hist_offset_ROSS.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))
                            if -180. <= np.mean(lon[step - 10:step + 10]) <= -130.:
                                #print('Ross')
                                offset_ROSS.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))
                                hist_offset_ROSS.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))
                            if -130. <= np.mean(lon[step - 10:step + 10]) <= -60.:
                                #print('Amundsen-Bellingshausen')
                                offset_AMBEL.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))
                                hist_offset_AMBEL.append(np.mean(ssha[step:step + 10]) - np.mean(ssha[step - 10:step]))

            monthly_offset_WEDD.append(np.mean(offset_WEDD))
            monthly_offset_IND.append(np.mean(offset_IND))
            monthly_offset_ROSS.append(np.mean(offset_ROSS))  
            monthly_offset_AMBEL.append(np.mean(offset_AMBEL))
            monthly_offset.append((np.mean(offset_WEDD) + np.mean(offset_IND)
                + np.mean(offset_ROSS) + np.mean(offset_AMBEL))/4)

    pl.figure()
    pl.hist([hist_offset_WEDD, hist_offset_IND, hist_offset_ROSS, hist_offset_AMBEL], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
    pl.legend()
    pl.title(yr + ' offset histogram')
    pl.xlabel('Offset Bin (m)')
    pl.ylabel('Frequency')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Quality Control/' + yr + '_offset_hist.png', format='png')
    pl.close()

    A = range(np.min(month_number), np.max(month_number) + 1)
    pl.figure()
    pl.plot(A, monthly_offset, label='Circumpolar')
    pl.plot(A, monthly_offset_WEDD, label='Weddell')
    pl.plot(A, monthly_offset_IND, label='Indian')
    pl.plot(A, monthly_offset_ROSS, label='Ross')
    pl.plot(A, monthly_offset_AMBEL, label='Amundsen-Bellingshausen')
    pl.legend(loc='lower right')
    pl.title(yr + ' monthly offset between lead and ocean ssha (m)')
    pl.ylabel('ssha offset (m)')
    pl.xlabel('Month')
    pl.xlim([1, 12])
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Quality Control/' + yr + '_monthly_offset_sectors.png', format='png')
    pl.close()