from exodus_read_base2 import ExodusRead
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit

out = ['01']
#out = ['01', '02']
num = 0

def func(x, a, b):
    return(a * np.exp(-x/b))


colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
num = 0

ymin = 0
ymax = 0
for v in out:
    exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_gopalakrishnan_ad_m1_kV_1um_'+v+'.e')

    plt_num = '001'


    cmap = plt.get_cmap('jet')

    xgrid = exopy.exodus.variables['coordx'][:]
    block0 = exopy.split_block('x', 0)
    block1 = exopy.split_block('x', 1)
    x0 = xgrid[block0]
    x1 = xgrid[block1]

    time_scale = 1

    dt_start = 0
    dt_end = -1

    gas_species = ['Arp', 'em']
    liquid_species = ['emliq', 'OHm', 'Om', 'O2m', 'O3m', 'HO2m', 'H+', 
                  'O', 'O2_1', 'O3', 'H', 'H2', 'HO2', 'OH', 'H2O2']

    #print(exopy.exodus.variables.keys())

    cell_block0 = exopy.exodus.variables['connect1'][:]
    cell_block1 = exopy.exodus.variables['connect2'][:]

    total_connect = np.vstack((cell_block0, cell_block1))
    num_cells = len(total_connect)


    phi = exopy.get_vals('potential', data_type='node')[:]

    efield_b0 = np.gradient(phi[-1,block0]*1000., x0)
    efield_b1 = np.gradient(phi[-1,block1]*1000., x1)

    plt.plot(x0, efield_b0, color=colors[num])
    plt.plot(x1, efield_b1, color=colors[num])

    ymax_temp = phi[-1, block0[-1]]
    ymin_temp = phi[-1, block1[0]]
    if (ymax_temp > ymax):
        ymax = np.copy(ymax_temp)
    if (ymin_temp < ymin):
        ymin = np.copy(ymin_temp)
    num += 1

#plt.xlim([1e-3 - 1e-6, 1e-3 + 1e-6])
#plt.ylim([ymin*1.05, ymax*1.05])
#plt.ylim([
plt.show()
