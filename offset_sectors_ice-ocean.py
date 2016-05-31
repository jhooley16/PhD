import os
import functions as funct
import numpy as np
import matplotlib.pyplot as pl

average_circumpolar_offset = np.zeros(12)

for yr in ['2011', '2012', '2013', '2014', '2015']:
    print(yr)
    
    monthly_offset_WEDD = []
    monthly_offset_IND = []
    monthly_offset_ROSS = []
    monthly_offset_AMBEL = []
    monthly_offset = []

    hist_offset_WEDD = []
    hist_offset_IND = []
    hist_offset_ROSS = []
    hist_offset_AMBEL = []
    hist_offset_circum = []

    month_number = []

    for mnth in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:

        offset_WEDD = []
        offset_IND = []
        offset_ROSS = []
        offset_AMBEL = []
        circum_offset = []

        if os.path.isdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + mnth + '_elev'):
            os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + yr + mnth + '_elev')
            month_number.append(int(mnth))
        
            for file in os.listdir():
                
                # 200 points equals the decorrelation scale (50 km) 
                ice_edge = [1] * 50 + [2] * 50
                ice_edge2 = [2] * 50 + [1] * 50

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
                    if surface_type[it:it + len(ice_edge)] == ice_edge:
                        iedge.append(it + len(ice_edge)//2)
                    if surface_type[it:it + len(ice_edge)] == ice_edge2:
                        iedge.append(it + len(ice_edge)//2)
                
                for step in iedge:
                    if surface_type[step - len(ice_edge)//2:step] == [1] * (len(ice_edge)//2):
                        # Ocean to Ice step
                        if abs(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2])) < .5:
                            # # Calculate circumpolar offset
                            circum_offset.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                            hist_offset_circum.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                            # Choose what the sector the offset point lies within
                            if -60. <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= 0.:
                                #print('Weddell')
                                offset_WEDD.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                                hist_offset_WEDD.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                            if 0 <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= 160.:
                                #print('Indian')
                                offset_IND.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                                hist_offset_IND.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                            if 160. <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= 180.:
                                #print('Ross')
                                offset_ROSS.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                                hist_offset_ROSS.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                            if -180. <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= -130.:
                                #print('Ross')
                                offset_ROSS.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                                hist_offset_ROSS.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                            if -130. <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= -60.:
                                #print('Amundsen-Bellingshausen')
                                offset_AMBEL.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))
                                hist_offset_AMBEL.append(np.mean(ssha[step - len(ice_edge)//2:step]) - np.mean(ssha[step:step + len(ice_edge)//2]))


                    if surface_type[step - len(ice_edge)//2:step] == [2] * (len(ice_edge)//2):
                        # Ice to Ocean step
                        if abs(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step])) < .5:
                            # Calculate circumpolar offset
                            circum_offset.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                            hist_offset_circum.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                            # Choose what the sector the offset point lies within
                            if -60. <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= 0.:
                                #print('Weddell')
                                offset_WEDD.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                                hist_offset_WEDD.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                            if 0 <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= 160.:
                                #print('Indian')
                                offset_IND.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                                hist_offset_IND.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                            if 160. <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= 180.:
                                #print('Ross')
                                offset_ROSS.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                                hist_offset_ROSS.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                            if -180. <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= -130.:
                                #print('Ross')
                                offset_ROSS.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                                hist_offset_ROSS.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                            if -130. <= np.mean(lon[step - len(ice_edge)//2:step + len(ice_edge)//2]) <= -60.:
                                #print('Amundsen-Bellingshausen')
                                offset_AMBEL.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))
                                hist_offset_AMBEL.append(np.mean(ssha[step:step + len(ice_edge)//2]) - np.mean(ssha[step - len(ice_edge)//2:step]))

            monthly_offset_WEDD.append(np.mean(offset_WEDD))
            monthly_offset_IND.append(np.mean(offset_IND))
            monthly_offset_ROSS.append(np.mean(offset_ROSS))  
            monthly_offset_AMBEL.append(np.mean(offset_AMBEL))
            monthly_offset.append(np.mean(circum_offset))

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
    
    average_circumpolar_offset += np.array(monthly_offset)

average_circumpolar_offset /= 5

print(average_circumpolar_offset)

print('January offset: ', average_circumpolar_offset[0])
print('Febuary offset: ', average_circumpolar_offset[1])
print('March offset: ', average_circumpolar_offset[2])
print('April offset: ', average_circumpolar_offset[3])
print('May offset: ', average_circumpolar_offset[4])
print('June offset: ', average_circumpolar_offset[5])
print('July offset: ', average_circumpolar_offset[6])
print('August offset: ', average_circumpolar_offset[7])
print('September offset: ', average_circumpolar_offset[8])
print('October offset: ', average_circumpolar_offset[9])
print('November offset: ', average_circumpolar_offset[10])
print('December offset: ', average_circumpolar_offset[11])

#January offset:  0.0532795957319
#Febuary offset:  0.0424895753339
#March offset:  0.044196774024
#April offset:  0.058855642473
#May offset:  0.0623327136638
#June offset:  0.0885369699514
#July offset:  0.108297693989
#August offset:  0.104956289149
#September offset:  0.0924955874413
#October offset:  0.0964708361887
#November offset:  0.0591056776193
#December offset:  0.0469972667101

