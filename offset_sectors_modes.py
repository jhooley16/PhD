import os
import functions as funct
import numpy as np
import matplotlib.pyplot as pl

average_circumpolar_offset_SAR = np.zeros(12)
average_circumpolar_offset_SARIn = np.zeros(12)

for year in ['2011', '2012', '2013', '2014', '2015']:
    print(year)
    
    monthly_offset_WEDD_SAR = []
    monthly_offset_IND_SAR = []
    monthly_offset_ROSS_SAR = []
    monthly_offset_AMBEL_SAR = []
    monthly_offset_SAR = []
    monthly_offset_WEDD_SARIn = []
    monthly_offset_IND_SARIn = []
    monthly_offset_ROSS_SARIn = []
    monthly_offset_AMBEL_SARIn = []
    monthly_offset_SARIn = []

    hist_offset_WEDD_SAR = []
    hist_offset_IND_SAR = []
    hist_offset_ROSS_SAR = []
    hist_offset_AMBEL_SAR = []
    hist_offset_circum_SAR = []
    hist_offset_WEDD_SARIn = []
    hist_offset_IND_SARIn = []
    hist_offset_ROSS_SARIn = []
    hist_offset_AMBEL_SARIn = []
    hist_offset_circum_SARIn = []

    month_number = []

    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:

        offset_WEDD_SAR = []
        offset_IND_SAR = []
        offset_ROSS_SAR = []
        offset_AMBEL_SAR = []
        circum_offset_SAR = []
        offset_WEDD_SARIn = []
        offset_IND_SARIn = []
        offset_ROSS_SARIn = []
        offset_AMBEL_SARIn = []
        circum_offset_SARIn = []

        if os.path.isdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev'):
            os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev')
            month_number.append(int(month))
        
            for file in os.listdir():
                
                # 200 points equals the decorrelation scale (50 km) 
                SAR_edge = [0] * 50 + [1] * 50
                SAR_edge2 = [1] * 50 + [0] * 50
                SARIn_edge = [1] * 50 + [2] * 50
                SARIn_edge2 = [2] * 50 + [1] * 50

                ssha = []
                lon = []
                lat = []
                surface_type = []

                f = open(file, 'r')
                for line in f:
                    line = line.strip()
                    columns = line.split()
                    # If data point is from open ocean (1) or from a lead (2)
                    if columns[0] == '1' or columns[0] == '2':
                        # If data point is listed as 'valid'
                        if columns[1] == '1':
                            lat.append(float(columns[5]))
                            lon.append(float(columns[6]))
                            ssha.append(float(columns[7]) - float(columns[8]))
                
                # Generate a list of retracker modes for this track
                surface_type = funct.mode_points(lat, lon, month)

                iedge_SAR = []
                iedge_SARIn = []
                for it in range(len(surface_type)):
                    if surface_type[it:it + len(SAR_edge)] == SAR_edge:
                        iedge_SAR.append(it + len(SAR_edge)//2)
                    elif surface_type[it:it + len(SAR_edge)] == SAR_edge2:
                        iedge_SAR.append(it + len(SAR_edge)//2)

                    elif surface_type[it:it + len(SARIn_edge)] == SARIn_edge:
                        iedge_SARIn.append(it + len(SARIn_edge)//2)
                    elif surface_type[it:it + len(SARIn_edge)] == SARIn_edge2:
                        iedge_SARIn.append(it + len(SARIn_edge)//2)

                for step in iedge_SAR:
                    if surface_type[step - len(SAR_edge)//2:step] == [0] * (len(SAR_edge)//2):
                        # LRM to SAR step
                        if abs(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2])) < .5:
                            # # Calculate circumpolar offset
                            circum_offset_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                            hist_offset_circum_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                            # Choose what the sector the offset point lies within
                            if -60. <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= 0.:
                                #print('Weddell')
                                offset_WEDD_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                                hist_offset_WEDD_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                            elif 0 <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= 160.:
                                #print('Indian')
                                offset_IND_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                                hist_offset_IND_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                            elif 160. <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= 180.:
                                #print('Ross')
                                offset_ROSS_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                                hist_offset_ROSS_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                            elif -180. <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= -130.:
                                #print('Ross')
                                offset_ROSS_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                                hist_offset_ROSS_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                            elif -130. <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= -60.:
                                #print('Amundsen-Bellingshausen')
                                offset_AMBEL_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))
                                hist_offset_AMBEL_SAR.append(np.mean(ssha[step - len(SAR_edge)//2:step]) - np.mean(ssha[step:step + len(SAR_edge)//2]))


                    if surface_type[step - len(SAR_edge)//2:step] == [1] * (len(SAR_edge)//2):
                        # SAR to LRM step
                        if abs(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step])) < .5:
                            # Calculate circumpolar offset
                            circum_offset_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                            hist_offset_circum_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                            # Choose what the sector the offset point lies within
                            if -60. <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= 0.:
                                #print('Weddell')
                                offset_WEDD_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                                hist_offset_WEDD_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                            elif 0 <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= 160.:
                                #print('Indian')
                                offset_IND_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                                hist_offset_IND_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                            elif 160. <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= 180.:
                                #print('Ross')
                                offset_ROSS_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                                hist_offset_ROSS_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                            elif -180. <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= -130.:
                                #print('Ross')
                                offset_ROSS_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                                hist_offset_ROSS_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                            elif -130. <= np.mean(lon[step - len(SAR_edge)//2:step + len(SAR_edge)//2]) <= -60.:
                                #print('Amundsen-Bellingshausen')
                                offset_AMBEL_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))
                                hist_offset_AMBEL_SAR.append(np.mean(ssha[step:step + len(SAR_edge)//2]) - np.mean(ssha[step - len(SAR_edge)//2:step]))

                for step in iedge_SARIn:
                    if surface_type[step - len(SAR_edge)//2:step] == [1] * (len(SAR_edge)//2):
                        # LRM to SAR step
                        if abs(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2])) < .5:
                            # # Calculate circumpolar offset
                            circum_offset_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                            hist_offset_circum_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                            # Choose what the sector the offset point lies within
                            if -60. <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= 0.:
                                #print('Weddell')
                                offset_WEDD_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                                hist_offset_WEDD_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                            elif 0 <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= 160.:
                                #print('Indian')
                                offset_IND_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                                hist_offset_IND_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                            elif 160. <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= 180.:
                                #print('Ross')
                                offset_ROSS_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                                hist_offset_ROSS_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                            elif -180. <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= -130.:
                                #print('Ross')
                                offset_ROSS_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                                hist_offset_ROSS_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                            elif -130. <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= -60.:
                                #print('Amundsen-Bellingshausen')
                                offset_AMBEL_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))
                                hist_offset_AMBEL_SARIn.append(np.mean(ssha[step - len(SARIn_edge)//2:step]) - np.mean(ssha[step:step + len(SARIn_edge)//2]))


                    if surface_type[step - len(SARIn_edge)//2:step] == [2] * (len(SARIn_edge)//2):
                        # SARIn to SAR step
                        if abs(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step])) < .5:
                            # Calculate circumpolar offset
                            circum_offset_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                            hist_offset_circum_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                            # Choose what the sector the offset point lies within
                            if -60. <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= 0.:
                                #print('Weddell')
                                offset_WEDD_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                                hist_offset_WEDD_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                            elif 0 <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= 160.:
                                #print('Indian')
                                offset_IND_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                                hist_offset_IND_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                            elif 160. <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= 180.:
                                #print('Ross')
                                offset_ROSS_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                                hist_offset_ROSS_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                            elif -180. <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= -130.:
                                #print('Ross')
                                offset_ROSS_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                                hist_offset_ROSS_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                            elif -130. <= np.mean(lon[step - len(SARIn_edge)//2:step + len(SARIn_edge)//2]) <= -60.:
                                #print('Amundsen-Bellingshausen')
                                offset_AMBEL_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))
                                hist_offset_AMBEL_SARIn.append(np.mean(ssha[step:step + len(SARIn_edge)//2]) - np.mean(ssha[step - len(SARIn_edge)//2:step]))

            monthly_offset_WEDD_SAR.append(np.mean(offset_WEDD_SAR))
            monthly_offset_IND_SAR.append(np.mean(offset_IND_SAR))
            monthly_offset_ROSS_SAR.append(np.mean(offset_ROSS_SAR))  
            monthly_offset_AMBEL_SAR.append(np.mean(offset_AMBEL_SAR))
            monthly_offset_SAR.append(np.mean(circum_offset_SAR))

            monthly_offset_WEDD_SARIn.append(np.mean(offset_WEDD_SARIn))
            monthly_offset_IND_SARIn.append(np.mean(offset_IND_SARIn))
            monthly_offset_ROSS_SARIn.append(np.mean(offset_ROSS_SARIn))  
            monthly_offset_AMBEL_SARIn.append(np.mean(offset_AMBEL_SARIn))
            monthly_offset_SARIn.append(np.mean(circum_offset_SARIn))

    pl.figure()
    pl.hist([hist_offset_WEDD_SAR, hist_offset_IND_SAR, hist_offset_ROSS_SAR, hist_offset_AMBEL_SAR], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
    pl.legend()
    pl.title(year + ' LRM - SAR histogram')
    pl.xlabel('Offset Bin (m)')
    pl.ylabel('Frequency')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Quality Control/' + yr + '_LRMSAR_offset_hist.png', format='png')
    pl.close()

    pl.figure()
    pl.hist([hist_offset_WEDD_SARIn, hist_offset_IND_SARIn, hist_offset_ROSS_SARIn, hist_offset_AMBEL_SARIn], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
    pl.legend()
    pl.title(year + ' SAR - SARIn histogram')
    pl.xlabel('Offset Bin (m)')
    pl.ylabel('Frequency')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Quality Control/' + year + '_SARSARIn_offset_hist.png', format='png')
    pl.close()

    A = range(np.min(month_number), np.max(month_number) + 1)
    pl.figure()
    pl.plot(A, monthly_offset_SAR, label='Circumpolar')
    pl.plot(A, monthly_offset_WEDD_SAR, label='Weddell')
    pl.plot(A, monthly_offset_IND_SAR, label='Indian')
    pl.plot(A, monthly_offset_ROSS_SAR, label='Ross')
    pl.plot(A, monthly_offset_AMBEL_SAR, label='Amundsen-Bellingshausen')
    pl.legend(loc='lower right')
    pl.title(year + ' monthly offset LRM - SAR (m)')
    pl.ylabel('Offset (m)')
    pl.xlabel('Month')
    pl.xlim([1, 12])
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Quality Control/' + year + '_monthly_LRMSAR_offset_sectors.png', format='png')
    pl.close()
    
    pl.figure()
    pl.plot(A, monthly_offset_SARIn, label='Circumpolar')
    pl.plot(A, monthly_offset_WEDD_SARIn, label='Weddell')
    pl.plot(A, monthly_offset_IND_SARIn, label='Indian')
    pl.plot(A, monthly_offset_ROSS_SARIn, label='Ross')
    pl.plot(A, monthly_offset_AMBEL_SARIn, label='Amundsen-Bellingshausen')
    pl.legend(loc='lower right')
    pl.title(year + ' monthly offset SAR - SARIn (m)')
    pl.ylabel('Offset (m)')
    pl.xlabel('Month')
    pl.xlim([1, 12])
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Quality Control/' + year + '_monthly_SARSARIn_offset_sectors.png', format='png')
    pl.close()

    average_circumpolar_offset_SAR += np.array(monthly_offset_SAR)
    average_circumpolar_offset_SARIn += np.array(monthly_offset_SARIn)

average_circumpolar_offset_SAR /= 5
average_circumpolar_offset_SARIn /= 5

print(average_circumpolar_offset_SAR)
print(average_circumpolar_offset_SARIn)

print('January LRM-SAR offset: ', average_circumpolar_offset_SAR[0])
print('Febuary LRM-SAR offset: ', average_circumpolar_offset_SAR[1])
print('March LRM-SAR offset: ', average_circumpolar_offset_SAR[2])
print('April LRM-SAR offset: ', average_circumpolar_offset_SAR[3])
print('May LRM-SAR offset: ', average_circumpolar_offset_SAR[4])
print('June LRM-SAR offset: ', average_circumpolar_offset_SAR[5])
print('July LRM-SAR offset: ', average_circumpolar_offset_SAR[6])
print('August LRM-SAR offset: ', average_circumpolar_offset_SAR[7])
print('September LRM-SAR offset: ', average_circumpolar_offset_SAR[8])
print('October LRM-SAR offset: ', average_circumpolar_offset_SAR[9])
print('November LRM-SAR offset: ', average_circumpolar_offset_SAR[10])
print('December LRM-SAR offset: ', average_circumpolar_offset_SAR[11])

print('January SAR-SARIn offset: ', average_circumpolar_offset_SARIn[0])
print('Febuary SAR-SARIn offset: ', average_circumpolar_offset_SARIn[1])
print('March SAR-SARIn offset: ', average_circumpolar_offset_SARIn[2])
print('April SAR-SARIn offset: ', average_circumpolar_offset_SARIn[3])
print('May SAR-SARIn offset: ', average_circumpolar_offset_SARIn[4])
print('June SAR-SARIn offset: ', average_circumpolar_offset_SARIn[5])
print('July SAR-SARIn offset: ', average_circumpolar_offset_SARIn[6])
print('August SAR-SARIn offset: ', average_circumpolar_offset_SARIn[7])
print('September SAR-SARIn offset: ', average_circumpolar_offset_SARIn[8])
print('October SAR-SARIn offset: ', average_circumpolar_offset_SARIn[9])
print('November SAR-SARIn offset: ', average_circumpolar_offset_SARIn[10])
print('December SAR-SARIn offset: ', average_circumpolar_offset_SARIn[11])

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

