import os
import functions as funct
import numpy as np
import matplotlib.pyplot as pl

average_circumpolar_offset_LRM_SAR = np.zeros(12)

WEDD_timeseries_LRM_SAR = []
IND_timeseries_LRM_SAR = []
ROSS_timeseries_LRM_SAR = []
AMBEL_timeseries_LRM_SAR = []
timeseries_LRM_SAR = []


for year in ['2011', '2012', '2013', '2014', '2015']:
    print(year)
    
    monthly_offset_WEDD_LRM_SAR = []
    monthly_offset_IND_LRM_SAR = []
    monthly_offset_ROSS_LRM_SAR = []
    monthly_offset_AMBEL_LRM_SAR = []
    monthly_offset_LRM_SAR = []

    hist_offset_WEDD_LRM_SAR = []
    hist_offset_IND_LRM_SAR = []
    hist_offset_ROSS_LRM_SAR = []
    hist_offset_AMBEL_LRM_SAR = []
    hist_offset_circum_LRM_SAR = []

    month_number = []

    for month in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
        print(month)
        offset_WEDD_LRM_SAR = []
        offset_IND_LRM_SAR = []
        offset_ROSS_LRM_SAR = []
        offset_AMBEL_LRM_SAR = []
        circum_offset_LRM_SAR = []

        if os.path.isdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev'):
            os.chdir('/Users/jmh2g09/Documents/PhD/Data/elev_files/' + year + month + '_elev')
            month_number.append(int(month))
        
            for file in os.listdir():
                
                # 200 points equals the decorrelation scale (50 km)
                SC = 100
                LRM_to_SAR = [0] * SC + [1] * SC
                SAR_to_LRM = [1] * SC + [0] * SC
                
                len_bound = len(LRM_to_SAR)             
                
                ssha = []
                lon = []
                lat = []
                surface = []

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
                            surface.append(int(columns[0]))
                
                # Generate a list of retracker modes for this track
                tracker_type = funct.mode_points(lat, lon, month)

                # Find the boundaries
                iedge_LRM_SAR = []
                for it in range(len(tracker_type)):
                    # Find retracker boundaries
                    if tracker_type[it:it + len_bound] == LRM_to_SAR:
                        iedge_LRM_SAR.append(it + len_bound//2)
                    elif tracker_type[it:it + len_bound] == SAR_to_LRM:
                        iedge_LRM_SAR.append(it + len_bound//2)

                # For the LRM-SAR boundary
                for step in iedge_LRM_SAR:
                    # If the boundary is ALL OCEAN
                    if surface[step - len_bound//2:step + len_bound//2] == [1] * len_bound:
                        # If it's a LRM to SAR step
                        if tracker_type[step - len_bound//2:step + len_bound//2] == LRM_to_SAR:
                            if np.max(abs(np.gradient(lat[step - len_bound//2:step + len_bound//2]))) < 0.05:
                                if abs(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2])) < .3:
                                    # Calculate circumpolar offset
                                    circum_offset_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    hist_offset_circum_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    # Choose what the sector the offset point lies within
                                    if -60. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 0.:
                                        #print('Weddell')
                                        offset_WEDD_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_WEDD_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    elif 0 <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 160.:
                                        #print('Indian')
                                        offset_IND_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_IND_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    elif 160. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 180.:
                                        #print('Ross')
                                        offset_ROSS_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_ROSS_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    elif -180. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -130.:
                                        #print('Ross')
                                        offset_ROSS_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_ROSS_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                    elif -130. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -60.:
                                        #print('Amundsen-Bellingshausen')
                                        offset_AMBEL_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))
                                        hist_offset_AMBEL_LRM_SAR.append(np.mean(ssha[step - len_bound//2:step]) - np.mean(ssha[step:step + len_bound//2]))

                        # If it's a SAR to LRM step
                        if tracker_type[step - len_bound//2:step + len_bound//2] == SAR_to_LRM:
                            if np.max(abs(np.gradient(lat[step - len_bound//2:step + len_bound//2]))) < 0.05:
                                if abs(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step])) < .3:
                                    # Calculate circumpolar offset
                                    circum_offset_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    hist_offset_circum_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    # Choose what the sector the offset point lies within
                                    if -60. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 0.:
                                        #print('Weddell')
                                        offset_WEDD_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_WEDD_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    elif 0 <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 160.:
                                        #print('Indian')
                                        offset_IND_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_IND_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    elif 160. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= 180.:
                                        #print('Ross')
                                        offset_ROSS_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_ROSS_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    elif -180. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -130.:
                                        #print('Ross')
                                        offset_ROSS_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_ROSS_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                    elif -130. <= np.mean(lon[step - len_bound//2:step + len_bound//2]) <= -60.:
                                        #print('Amundsen-Bellingshausen')
                                        offset_AMBEL_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))
                                        hist_offset_AMBEL_LRM_SAR.append(np.mean(ssha[step:step + len_bound//2]) - np.mean(ssha[step - len_bound//2:step]))

            monthly_offset_WEDD_LRM_SAR.append(np.mean(offset_WEDD_LRM_SAR))
            monthly_offset_IND_LRM_SAR.append(np.mean(offset_IND_LRM_SAR))
            monthly_offset_ROSS_LRM_SAR.append(np.mean(offset_ROSS_LRM_SAR))  
            monthly_offset_AMBEL_LRM_SAR.append(np.mean(offset_AMBEL_LRM_SAR))
            monthly_offset_LRM_SAR.append(np.mean(circum_offset_LRM_SAR))

            WEDD_timeseries_LRM_SAR.append(np.mean(offset_WEDD_LRM_SAR))
            IND_timeseries_LRM_SAR.append(np.mean(offset_IND_LRM_SAR))
            ROSS_timeseries_LRM_SAR.append(np.mean(offset_ROSS_LRM_SAR)) 
            AMBEL_timeseries_LRM_SAR.append(np.mean(offset_AMBEL_LRM_SAR))
            timeseries_LRM_SAR.append(np.mean(circum_offset_LRM_SAR))
    
    pl.figure()
    pl.hist([hist_offset_WEDD_LRM_SAR, hist_offset_IND_LRM_SAR, hist_offset_ROSS_LRM_SAR, hist_offset_AMBEL_LRM_SAR], label=['Weddell', 'Indian', 'Ross', 'Amundsen-Bellingshausen'])
    pl.legend(loc='top left')
    pl.title(year + ' LRM - SAR histogram')
    pl.xlabel('Offset Bin (m)')
    pl.ylabel('Frequency')
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_LRM_SAR_offset_hist.png', format='png')
    pl.close()

    A = range(np.min(month_number), np.max(month_number) + 1)
    pl.figure()
    pl.plot(A, monthly_offset_LRM_SAR, label='Circumpolar', marker='.')
    pl.plot(A, monthly_offset_WEDD_LRM_SAR, label='Weddell', marker='.')
    pl.plot(A, monthly_offset_IND_LRM_SAR, label='Indian', marker='.')
    pl.plot(A, monthly_offset_ROSS_LRM_SAR, label='Ross', marker='.')
    pl.plot(A, monthly_offset_AMBEL_LRM_SAR, label='Amundsen-Bellingshausen', marker='.')
    pl.legend(loc='top')
    pl.title(year + ' LRM - SAR (m) offset')
    pl.ylabel('Offset (m)')
    pl.xlabel('Month')
    pl.xlim([1, 12])
    pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/' + year + '_LRM_SAR_offset_sectors.png', format='png')
    pl.close()
    
    average_circumpolar_offset_LRM_SAR += np.array(monthly_offset_LRM_SAR)

f=open('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/LRM_SAR_timeseries.txt', 'w')
print(WEDD_timeseries_LRM_SAR, IND_timeseries_LRM_SAR, ROSS_timeseries_LRM_SAR, AMBEL_timeseries_LRM_SAR, timeseries_LRM_SAR, file=f)
f.close()

pl.figure()
pl.plot(timeseries_LRM_SAR, label='Circumpolar', marker='.')
pl.plot(WEDD_timeseries_LRM_SAR, label='Weddell', marker='.')
pl.plot(IND_timeseries_LRM_SAR, label='Indian', marker='.')
pl.plot(ROSS_timeseries_LRM_SAR, label='Ross', marker='.')
pl.plot(AMBEL_timeseries_LRM_SAR, label='Amundsen-Bellingshausen', marker='.')
pl.title('LRM - SAR (m) offset timeseries')
pl.xlabel('Month (from Jan 2011)')
pl.ylabel('LRM - SAR (m) offset')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/Seperate Modes/Figures/LRMSAR_timeseries.png', format='png')
pl.close()

average_circumpolar_offset_LRM_SAR /= 5

print(average_circumpolar_offset_LRM_SAR)

print('January LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[0])
print('Febuary LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[1])
print('March LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[2])
print('April LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[3])
print('May LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[4])
print('June LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[5])
print('July LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[6])
print('August LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[7])
print('September LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[8])
print('October LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[9])
print('November LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[10])
print('December LRM-SAR offset: ', average_circumpolar_offset_LRM_SAR[11])

#January LRM-SAR offset:  0.00489492491322
#Febuary LRM-SAR offset:  0.00688199790682
#March LRM-SAR offset:  -0.00897106525092
#April LRM-SAR offset:  -0.0181416179892
#May LRM-SAR offset:  -0.0101094678668
#June LRM-SAR offset:  -0.0187603584627
#July LRM-SAR offset:  -0.0295892239469
#August LRM-SAR offset:  -0.0219179266994
#September LRM-SAR offset:  -0.0162279361396
#October LRM-SAR offset:  -0.0174250239166
#November LRM-SAR offset:  0.00133077159719
#December LRM-SAR offset:  0.0124794197335
