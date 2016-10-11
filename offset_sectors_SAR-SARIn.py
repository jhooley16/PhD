import os
import functions as funct
import numpy as np
import matplotlib.pyplot as pl

average_circumpolar_offset_SAR = np.zeros(12)
average_circumpolar_offset_SARIn = np.zeros(12)

WEDD_timeseries_SAR = []
IND_timeseries_SAR = []
ROSS_timeseries_SAR = []
AMBEL_timeseries_SAR = []
timeseries_SAR = []
WEDD_timeseries_SARIn = []
IND_timeseries_SARIn = []
ROSS_timeseries_SARIn = []
AMBEL_timeseries_SARIn = []
timeseries_SARIn = []


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
        print(month)
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
                SAR_edge = [0] * 200 + [1] * 200
                SAR_edge2 = [1] * 200 + [0] * 200
                SARIn_edge = [1] * 200 + [2] * 200
                SARIn_edge2 = [2] * 200 + [1] * 200
                Len_SAR = len(SAR_edge)
                Len_SARIn = len(SARIn_edge)        
                
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
                    if surface_type[it:it + Len_SAR] == SAR_edge:
                        iedge_SAR.append(it + Len_SAR//2)
                    elif surface_type[it:it + Len_SAR] == SAR_edge2:
                        iedge_SAR.append(it + Len_SAR//2)

                    elif surface_type[it:it + Len_SARIn] == SARIn_edge:
                        iedge_SARIn.append(it + Len_SARIn//2)
                    elif surface_type[it:it + Len_SARIn] == SARIn_edge2:
                        iedge_SARIn.append(it + Len_SARIn//2)

                for step in iedge_SAR:
                    if surface_type[step - Len_SAR//2:step] == [0] * (Len_SAR//2):
                        # LRM to SAR step
                        if np.max(abs(np.gradient(lat[step - Len_SAR//2:step + Len_SAR//2]))) < 0.05:
                            if abs(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2])) < .3:
                                # # Calculate circumpolar offset
                                circum_offset_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                hist_offset_circum_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                # Choose what the sector the offset point lies within
                                if -60. <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= 0.:
                                    #print('Weddell')
                                    offset_WEDD_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                    hist_offset_WEDD_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                elif 0 <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= 160.:
                                    #print('Indian')
                                    offset_IND_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                    hist_offset_IND_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                elif 160. <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= 180.:
                                    #print('Ross')
                                    offset_ROSS_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                    hist_offset_ROSS_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                elif -180. <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= -130.:
                                    #print('Ross')
                                    offset_ROSS_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                    hist_offset_ROSS_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                elif -130. <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= -60.:
                                    #print('Amundsen-Bellingshausen')
                                    offset_AMBEL_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))
                                    hist_offset_AMBEL_SAR.append(np.mean(ssha[step - Len_SAR//2:step]) - np.mean(ssha[step:step + Len_SAR//2]))


                    if surface_type[step - Len_SAR//2:step] == [1] * (Len_SAR//2):
                        # SAR to LRM step
                        if np.max(abs(np.gradient(lat[step - Len_SAR//2:step + Len_SAR//2]))) < 0.05:
                            if abs(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step])) < .3:
                                # Calculate circumpolar offset
                                circum_offset_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                hist_offset_circum_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                # Choose what the sector the offset point lies within
                                if -60. <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= 0.:
                                    #print('Weddell')
                                    offset_WEDD_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                    hist_offset_WEDD_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                elif 0 <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= 160.:
                                    #print('Indian')
                                    offset_IND_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                    hist_offset_IND_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                elif 160. <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= 180.:
                                    #print('Ross')
                                    offset_ROSS_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                    hist_offset_ROSS_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                elif -180. <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= -130.:
                                    #print('Ross')
                                    offset_ROSS_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                    hist_offset_ROSS_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                elif -130. <= np.mean(lon[step - Len_SAR//2:step + Len_SAR//2]) <= -60.:
                                    #print('Amundsen-Bellingshausen')
                                    offset_AMBEL_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))
                                    hist_offset_AMBEL_SAR.append(np.mean(ssha[step:step + Len_SAR//2]) - np.mean(ssha[step - Len_SAR//2:step]))

                for step in iedge_SARIn:
                    if surface_type[step - Len_SARIn//2:step] == [1] * (Len_SARIn//2):
                        if np.max(abs(np.gradient(lat[step - Len_SARIn//2:step + Len_SARIn//2]))) < 0.05:
                            # SAR to SARIn step
                            if abs(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2])) < .3:
                                # # Calculate circumpolar offset
                                circum_offset_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                hist_offset_circum_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                # Choose what the sector the offset point lies within
                                if -60. <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= 0.:
                                    #print('Weddell')
                                    offset_WEDD_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                    hist_offset_WEDD_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                elif 0 <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= 160.:
                                    #print('Indian')
                                    offset_IND_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                    hist_offset_IND_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                elif 160. <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= 180.:
                                    #print('Ross')
                                    offset_ROSS_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                    hist_offset_ROSS_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                elif -180. <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= -130.:
                                    #print('Ross')
                                    offset_ROSS_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                    hist_offset_ROSS_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                elif -130. <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= -60.:
                                    #print('Amundsen-Bellingshausen')
                                    offset_AMBEL_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))
                                    hist_offset_AMBEL_SARIn.append(np.mean(ssha[step - Len_SARIn//2:step]) - np.mean(ssha[step:step + Len_SARIn//2]))


                    if surface_type[step - Len_SARIn//2:step] == [2] * (Len_SARIn//2):
                        if np.max(abs(np.gradient(lat[step - Len_SARIn//2:step + Len_SARIn//2]))) < 0.05:
                            # SARIn to SAR step
                            if abs(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step])) < .3:
                                # Calculate circumpolar offset
                                circum_offset_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                hist_offset_circum_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                # Choose what the sector the offset point lies within
                                if -60. <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= 0.:
                                    #print('Weddell')
                                    offset_WEDD_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                    hist_offset_WEDD_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                elif 0 <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= 160.:
                                    #print('Indian')
                                    offset_IND_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                    hist_offset_IND_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                elif 160. <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= 180.:
                                    #print('Ross')
                                    offset_ROSS_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                    hist_offset_ROSS_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                elif -180. <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= -130.:
                                    #print('Ross')
                                    offset_ROSS_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                    hist_offset_ROSS_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                elif -130. <= np.mean(lon[step - Len_SARIn//2:step + Len_SARIn//2]) <= -60.:
                                    #print('Amundsen-Bellingshausen')
                                    offset_AMBEL_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))
                                    hist_offset_AMBEL_SARIn.append(np.mean(ssha[step:step + Len_SARIn//2]) - np.mean(ssha[step - Len_SARIn//2:step]))

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

            WEDD_timeseries_SAR.append(np.mean(offset_WEDD_SAR))
            IND_timeseries_SAR.append(np.mean(offset_IND_SAR))
            ROSS_timeseries_SAR.append(np.mean(offset_ROSS_SAR)) 
            AMBEL_timeseries_SAR.append(np.mean(offset_AMBEL_SAR))
            timeseries_SAR.append(np.mean(circum_offset_SAR))
            WEDD_timeseries_SARIn.append(np.mean(offset_WEDD_SARIn))
            IND_timeseries_SARIn.append(np.mean(offset_IND_SARIn))
            ROSS_timeseries_SARIn.append(np.mean(offset_ROSS_SARIn))  
            AMBEL_timeseries_SARIn.append(np.mean(offset_AMBEL_SARIn))
            timeseries_SARIn.append(np.mean(circum_offset_SARIn))
    
    pl.figure()
    pl.hist([hist_offset_WEDD_SAR, hist_offset_IND_SAR, hist_offset_ROSS_SAR, hist_offset_AMBEL_SAR], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
    pl.legend()
    pl.title(year + ' LRM - SAR histogram')
    pl.xlabel('Offset Bin (m)')
    pl.ylabel('Frequency')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_LRMSAR_offset_hist.png', format='png')
    pl.close()

    pl.figure()
    pl.hist([hist_offset_WEDD_SARIn, hist_offset_IND_SARIn, hist_offset_ROSS_SARIn, hist_offset_AMBEL_SARIn], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
    pl.legend()
    pl.title(year + ' SAR - SARIn histogram')
    pl.xlabel('Offset Bin (m)')
    pl.ylabel('Frequency')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_SARSARIn_offset_hist.png', format='png')
    pl.close()

    A = range(np.min(month_number), np.max(month_number) + 1)
    pl.figure()
    pl.plot(A, monthly_offset_SAR, label='Circumpolar', marker='.')
    pl.plot(A, monthly_offset_WEDD_SAR, label='Weddell', marker='.')
    pl.plot(A, monthly_offset_IND_SAR, label='Indian', marker='.')
    pl.plot(A, monthly_offset_ROSS_SAR, label='Ross', marker='.')
    pl.plot(A, monthly_offset_AMBEL_SAR, label='Amundsen-Bellingshausen', marker='.')
    pl.legend(loc='lower right')
    pl.title(year + ' monthly offset LRM - SAR (m)')
    pl.ylabel('Offset (m)')
    pl.xlabel('Month')
    pl.xlim([1, 12])
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_monthly_LRMSAR_offset_sectors.png', format='png')
    pl.close()
    
    pl.figure()
    pl.plot(A, monthly_offset_SARIn, label='Circumpolar', marker='.')
    pl.plot(A, monthly_offset_WEDD_SARIn, label='Weddell', marker='.')
    pl.plot(A, monthly_offset_IND_SARIn, label='Indian', marker='.')
    pl.plot(A, monthly_offset_ROSS_SARIn, label='Ross', marker='.')
    pl.plot(A, monthly_offset_AMBEL_SARIn, label='Amundsen-Bellingshausen', marker='.')
    pl.legend(loc='lower right')
    pl.title(year + ' monthly offset SAR - SARIn (m)')
    pl.ylabel('Offset (m)')
    pl.xlabel('Month')
    pl.xlim([1, 12])
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_monthly_SARSARIn_offset_sectors.png', format='png')
    pl.close()
    
    average_circumpolar_offset_SAR += np.array(monthly_offset_SAR)
    average_circumpolar_offset_SARIn += np.array(monthly_offset_SARIn)

f=open('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/SAR_timeseries.txt', 'w')
for i in range(len(WEDD_timeseries_SAR)):
    print(WEDD_timeseries_SAR[i], IND_timeseries_SAR[i], ROSS_timeseries_SAR[i], AMBEL_timeseries_SAR[i], timeseries_SAR[i], file=f)
f.close()
f=open('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/SARIn_timeseries.txt', 'w')
for i in range(len(WEDD_timeseries_SARIn)):
    print(WEDD_timeseries_SARIn[i], IND_timeseries_SARIn[i], ROSS_timeseries_SARIn[i], AMBEL_timeseries_SARIn[i], timeseries_SARIn[i], file=f)
f.close()

pl.figure()
pl.plot(timeseries_SAR, label='Circumpolar', marker='.')
pl.plot(WEDD_timeseries_SAR, label='Weddell', marker='.')
pl.plot(IND_timeseries_SAR, label='Indian', marker='.')
pl.plot(ROSS_timeseries_SAR, label='Ross', marker='.')
pl.plot(AMBEL_timeseries_SAR, label='Amundsen-Bellingshausen', marker='.')
pl.title('LRM - SAR (m) offset timeseries')
pl.xlabel('Month (from Jan 2011)')
pl.ylabel('LRM - SAR (m) offset')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/LRMSAR_timeseries.png', format='png')
pl.close()

pl.figure()
pl.plot(timeseries_SARIn, label='Circumpolar', marker='.')
pl.plot(WEDD_timeseries_SARIn, label='Weddell', marker='.')
pl.plot(IND_timeseries_SARIn, label='Indian', marker='.')
pl.plot(ROSS_timeseries_SARIn, label='Ross', marker='.')
pl.plot(AMBEL_timeseries_SARIn, label='Amundsen-Bellingshausen', marker='.')
pl.title('SAR - SARIn (m) offset timeseries', marker='.')
pl.xlabel('Month (from Jan 2011)')
pl.ylabel('SAR - SARIn (m) offset')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/SARSARIn_timeseries.png', format='png')
pl.close()

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

#January LRM-SAR offset:  0.00605720510952
#Febuary LRM-SAR offset:  0.0077402495848
#March LRM-SAR offset:  -0.00528115908118
#April LRM-SAR offset:  -0.0122375772566
#May LRM-SAR offset:  -0.004982935448
#June LRM-SAR offset:  -0.0157234755177
#July LRM-SAR offset:  -0.0179909306721
#August LRM-SAR offset:  -0.0129872770328
#September LRM-SAR offset:  -0.0123001013372
#October LRM-SAR offset:  -0.0102008805335
#November LRM-SAR offset:  0.00215810350424
#December LRM-SAR offset:  0.0109285059187
#January SAR-SARIn offset:  -0.0010534507906
#Febuary SAR-SARIn offset:  -0.00611335649884
#March SAR-SARIn offset:  -0.00755736513066
#April SAR-SARIn offset:  -0.00960064774106
#May SAR-SARIn offset:  -0.00895712247207
#June SAR-SARIn offset:  -0.00773984717228
#July SAR-SARIn offset:  -0.00492843462298
#August SAR-SARIn offset:  -0.00431128796423
#September SAR-SARIn offset:  -0.00282061906754
#October SAR-SARIn offset:  -0.00117488669705
#November SAR-SARIn offset:  -0.00406919339594
#December SAR-SARIn offset:  -0.0011911462509

