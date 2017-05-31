import os
import functions as funct
import numpy as np
import matplotlib.pyplot as pl
from datetime import date

LRM_SAR_offset = []
ocean_ice_offset = []
for imnth in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
    LRM_SAR_offset = np.append(LRM_SAR_offset, funct.apply_offset(imnth, 'SAR_ocean'))
    ocean_ice_offset = np.append(ocean_ice_offset, funct.apply_offset(imnth, 'ice'))

LRM_SAR_constant = np.full(len(LRM_SAR_offset), fill_value=0) + funct.apply_offset('constant', 'SAR_ocean')
ocean_ice_constant = np.full(len(LRM_SAR_offset), fill_value=0) + funct.apply_offset('constant', 'ice')

pl.figure()
pl.grid()
pl.plot(range(1, 13), LRM_SAR_offset * 100, label='SAR$_{ocean}$', color='b')
pl.plot(range(1, 13), LRM_SAR_constant * 100, label='SAR$_{ocean}$ Constant $=$ ' + str(np.round(funct.apply_offset('constant', 'SAR_ocean') * 100, 1)) + ' cm', color='b', ls='--')
pl.plot(range(1, 13), ocean_ice_offset * 100, label='SAR$_{lead}$', color='r')
pl.plot(range(1, 13), ocean_ice_constant * 100, label='SAR$_{lead}$ Constant $=$ ' + str(np.round(funct.apply_offset('constant', 'ice') * 100, 1)) + ' cm', color='r', ls='--')
my_xticks = ['Jan','Feb','Mar','Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec']
pl.xticks(range(1, 13), my_xticks)
pl.xlim(1, 12)
pl.ylim(-5, 9)
pl.legend(loc='best', prop={'size':9})
pl.ylabel('Offset (cm)')
pl.savefig('/Users/jmh2g09/Documents/PhD/Data/SeparateModes/Figures/OFFSET.png', doi=300, transparent=True, bbox_inches='tight')
pl.close()