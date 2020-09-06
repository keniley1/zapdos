from exodus_read_base2 import ExodusRead
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import gridspec
from scipy.ndimage import gaussian_filter1d

def gaussian_function(x, y, sigma):
    '''
    x: grid
    y: current location on grid
    '''
    f = 1.0/np.sqrt(2.0*np.pi)/sigma * np.exp(-0.5*((x-y)/sigma)**2.0)
    return(f / np.trapz(f, x=x))

def average_over_width(x, y, window):
    '''
    Averages values of y(x) onto a coarse grid of width `window`.
    '''
    x_coarse = np.arange(x[0], x[-1], window)
    y_coarse = np.zeros(shape=(len(y), len(x_coarse)))
    if (x_coarse[-1] < x[-1]):
        x_temp = np.append(x_coarse, x[-1])
    else:
        x_temp = np.copy(x_coarse)

    for t in range(len(x_temp)-1):
        test = np.where(np.logical_and(x > x_temp[t], x < x_temp[t+1]))[:][0]
        #y_coarse[:,t] = np.mean(y[:,test])
        y_coarse[:,t] = np.trapz(y[:, test], x=x[test])/(x[test[-1]] - x[test[0]])
    return(x_coarse, y_coarse)

def step_function(x, y, sigma):
    '''
    Defines a step function. Essentially a equal-weight filter of finite width.
    '''
    f = 1.0/x
    return(f / np.trapz(f, x=x))

def gaussian_filter(f, x, sigma, ti):
    '''
    Accepts f (function) and x (grid)  as input.
    Applies gaussian filter (Weierstrass transform)
    Returns filtered data, F
    Only works on 1D data!
    '''
    F = np.zeros(shape=(np.shape(f)))
    for i in range(len(x)):
        F[ti,i] = np.trapz(gaussian_function(x, x[i], sigma) * f[ti,:], x=x) 
    return(F)

def time_range(time, t_range):
    if (len(t_range) != 2):
        print("ERROR: t_range should be a list with two values: the lower and upper bounds of time.")
        exit()
    return(np.where((time > t_range[0]) & (time < t_range[1]))[0][:])

def window_average(x, y, time_window):
    return(np.trapz(y, x=time_window)/(time_window[-1] - time_window[0]))

def line_average(x, y, time, voltage, time_ranges, voltage_window, species='', use_mole=True, mult_factor=1.0):
    fig = plt.figure(figsize=(10,8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1,3])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    tmin = min(min(time_ranges))
    tmax = max(max(time_ranges))

    if (use_mole):
        new_y = np.exp(y) * 6.022e23
    else:
        new_y = np.copy(y)

    # Compute average over time slice
    #yavg = np.sum(y[time_index,:], axis=0)/float(len(time_index))
    # (also convert to cm^-3)
    yavg = np.zeros(shape=(len(x),))
    for i in range(len(time_ranges)):
        time_index = time_range(time, time_ranges[i])
        val = np.swapaxes(new_y[time_index,:], 0, 1)
        yavg += np.trapz(val, x=time[time_index])/(time[time_index[-1]] - time[time_index[0]])
    yavg = yavg / float(len(time_ranges)) * mult_factor

    # Plot voltage (top plot)
    ax0.plot(time[voltage_window], voltage[voltage_window], 'k', linewidth=2)
    ax0.axvspan(tmin, tmax, alpha=0.5, color='red')


    ax1.plot(x, yavg, 'k', linewidth=2, label=species)

    return(fig,ax0,ax1)

#exopy = ExodusRead('/home/shane/projects/zapdos/problems/air_dry_gold_files/sim_2ms/air_dry_reduced_out_002.e')
exopy = ExodusRead('/home/shane/projects/zapdos/problems/air_dry_gold_files/sim_2ms/air_dry_added_exodus_out.e')


plt_num = '002'

validation_plot_loc = 'plt2/'
#validation_plot_loc = 'validation_plots/'


######################################################
# plt_num = '001'  -- air_dry_charge_vib_out.e
# plt_num = '002'  -- air_dry_charge_vib_out_RATES.e
######################################################

xgrid = exopy.exodus.variables['coordx'][:]

time_scale = 1e6

dt_start = 700
#dt_end = 3379
dt_end = 1500

species = ['O3', 'O21S', 'O2-', 'O2+', 'O2*', 'O-', 'O+', 'O*', 'O', 'NO3', 
           'NO2', 'NO+', 'NO', 'N4+', 'N2v', 'N2O5', 'N2O4', 'N2O3', 'N2+', 
           'N2**', 'N2*', 'N+', 'N*', 'N', 'em']

#print(exopy.exodus.variables.keys())

cell_block0 = exopy.exodus.variables['connect1'][:]
cell_block1 = exopy.exodus.variables['connect2'][:]
cell_block2 = exopy.exodus.variables['connect3'][:]

total_connect = np.vstack((cell_block0, cell_block1, cell_block2))
num_cells = len(total_connect)

block0 = exopy.split_block('x', 0)
block1 = exopy.split_block('x', 1)
block2 = exopy.split_block('x', 2)

xcell = np.zeros(shape=(num_cells))
for i in range(num_cells):
    xcell[i] = (xgrid[total_connect[i][1]-1] + xgrid[total_connect[i][0]-1])/2.0

xcell1 = np.zeros(shape=(len(cell_block1)))
for i in range(len(cell_block1)):
    xcell1[i] = (xgrid[cell_block1[i][1]-1] + xgrid[cell_block1[i][0]-1])/2.0

time = exopy.time[:]
# Note that time is in units of MICROSECONDS
first_three_halfwaves = [[0.7, 1.7], [1.7, 2.7], [2.7, 3.7]]
#first_three_halfwaves = [[0.7e-6, 3.7e-6]]
first_negative = [[0.7, 1.7]]
full_pulse = [[0.4, 5]]

#first_three_halfwaves = [[0.34, 1.34], [1.34, 2.34], [2.34, 3.34]]
#first_negative = [[0.34, 1.34]]
#first_three_halfwaves = [[1.34, 2.34], [2.34, 3.34], [3.34, 4.34]]
#first_negative = [[1.34, 2.34]]
first_three_halfwaves = [[0, 1.34], [1.34, 2.34], [2.34, 3.34]]
first_negative = [[0, 1.34]]
full_pulse = [[0.4, 5]]
#full_pulse = [[1.34, 5]]

tmin = 1e-1
#tmin = 1e-6
#tmax = 1.4e-6
#tmin = 0.6e-6
tmax = 3.7
# time_slice is the variable that dictates where the average will be taken place
time_slice = np.where((time > tmin) & (time < tmax))[0][:]
time_index = np.arange(len(exopy.time[time_slice]))

phi0 = exopy.get_vals('potential_dom0', data_type='node')[:]


voltage_window = np.where(time < 7e-6)[0][:]
Vw = phi0[voltage_window,0]

em = exopy.get_vals('em', data_type='node')[:]
phi = exopy.get_vals('potential_dom1', data_type='node')[:]

time_index = time_range(time*1e6, first_negative[0])
ftest = gaussian_filter(em[:, block1], xgrid[block1], 8.3e-6, time_index)
#fig1, ax0, ax1 = line_average(xgrid[block1], ftest, time*1e6, Vw, first_three_halfwaves, voltage_window, species='e$^-$', mult_factor=1e-6)
fig1, ax0, ax1 = line_average(xgrid[block1], ftest, time*1e6, Vw, first_negative, voltage_window, species='e$^-$', mult_factor=1e-6)
ax0.set_xlabel('Time ($\mu$s)', fontsize=16)
ax0.set_ylabel('Voltage [kV]', fontsize=16)
ax1.set_xlabel('Gap Distance (m)', fontsize=16)
ax1.set_ylabel('Density (cm$^{-3}$)', fontsize=16)
#ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax0.xaxis.set_label_position('top')
ax0.xaxis.tick_top()
ax0.tick_params(axis='both', labelsize=12)
ax1.tick_params(axis='both', labelsize=12)
ax1.legend(loc='upper right', frameon=False, fontsize=14)
plt.savefig(validation_plot_loc+'ne_first_halfwave_'+plt_num+'.png', dpi=200, bbox_inches='tight')
plt.close()


ftest = gaussian_filter(em[:, block1], xgrid[block1], 45e-6, time_index)
#fig1, ax0, ax1 = line_average(xgrid[block1], ftest, time*1e6, Vw, first_three_halfwaves, voltage_window, species='e$^-$', mult_factor=1e-6)
x1,y1 = average_over_width(xgrid[block1], em[:, block1], 45e-6)
fig1, ax0, ax1 = line_average(x1, y1, time*1e6, Vw, first_negative, voltage_window, species='e$^-$', mult_factor=1e-6)
ax0.set_xlabel('Time ($\mu$s)', fontsize=16)
ax0.set_ylabel('Voltage [kV]', fontsize=16)
ax1.set_xlabel('Gap Distance (m)', fontsize=16)
ax1.set_ylabel('Density (cm$^{-3}$)', fontsize=16)
#ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax0.xaxis.set_label_position('top')
ax0.xaxis.tick_top()
ax0.tick_params(axis='both', labelsize=12)
ax1.tick_params(axis='both', labelsize=12)
ax1.legend(loc='upper right', frameon=False, fontsize=14)
plt.close()
#plt.savefig('validation_plots/ne_first_three_halfwaves.png', bbox_inches='tight', dpi=200)
#plt.savefig('validation_plots/ne_first_negative.png', bbox_inches='tight', dpi=200)

# Electric field test
efield = exopy.get_vals('Efield', data_type='element', block=1)[:]
x1, y1 = average_over_width(xcell1, abs(efield)/2.245e25 * 1e21, 45e-6) 
fig2, ax20, ax21 = line_average(x1, y1, time*1e6, Vw, first_three_halfwaves, voltage_window, species='Reduced Field', use_mole=False) 
ax20.set_xlabel('Time ($\mu$s)', fontsize=16)
ax20.set_ylabel('Voltage [kV]', fontsize=16)
ax21.set_xlabel('Gap Distance (m)', fontsize=16)
ax21.set_ylabel('Reduced Field (Td)', fontsize=16)
ax20.xaxis.set_label_position('top')
ax20.xaxis.tick_top()
ax20.tick_params(axis='both', labelsize=12)
ax21.tick_params(axis='both', labelsize=12)
ax21.legend(loc='upper right', frameon=False, fontsize=14)
#plt.savefig('validation_plots/efield_first_three_halfwaves'+plt_num+'.png', bbox_inches='tight', dpi=200)
#plt.savefig('validation_plots/efield_first_negative'+plt_num+'.png', bbox_inches='tight', dpi=200)
plt.savefig(validation_plot_loc+'efield_first_halfwave_'+plt_num+'.png', bbox_inches='tight', dpi=200)
plt.close()

# ELECTRIC FIELD
efield = exopy.get_vals('Efield', data_type='element', block=1)[:]
fig1, ax0, ax1 = line_average(xgrid[block1], em[:, block1], time*1e6, Vw, first_three_halfwaves, voltage_window, species='e$^-$', mult_factor=1e-6)
#fig1, ax0, ax1 = line_average(xgrid[block1], em[:, block1], time*1e6, Vw, first_negative, voltage_window, species='e$^-$')
ax0.set_xlabel('Time ($\mu$s)', fontsize=16)
ax0.set_ylabel('Voltage [kV]', fontsize=16)
ax1.set_xlabel('Gap Distance (m)', fontsize=16)
ax1.set_ylabel('Density (cm$^{-3}$)', fontsize=16)
#ax0.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax0.xaxis.set_label_position('top')
ax0.xaxis.tick_top()
ax0.tick_params(axis='both', labelsize=12)
ax1.tick_params(axis='both', labelsize=12)
ax1.legend(loc='upper right', frameon=False, fontsize=14)
plt.savefig(validation_plot_loc+'ne_first_three_halfwaves_'+plt_num+'.png', bbox_inches='tight', dpi=200)
#plt.savefig('validation_plots/ne_first_negative_'+plt_num+'.png', bbox_inches='tight', dpi=200)
plt.close()


fig2, ax20, ax21 = line_average(xcell1, abs(efield)/2.245e25*1e21, time*1e6, Vw, first_three_halfwaves, voltage_window, species='Reduced Field', use_mole=False) 
#fig2, ax20, ax21 = line_average(xcell1, abs(efield)/2.245e25*1e21, time*1e6, Vw, first_negative, voltage_window, species='Reduced Field', use_mole=False) 
ax20.set_xlabel('Time ($\mu$s)', fontsize=16)
ax20.set_ylabel('Voltage [kV]', fontsize=16)
ax21.set_xlabel('Gap Distance (m)', fontsize=16)
ax21.set_ylabel('Reduced Field (Td)', fontsize=16)
ax20.xaxis.set_label_position('top')
ax20.xaxis.tick_top()
ax20.tick_params(axis='both', labelsize=12)
ax21.tick_params(axis='both', labelsize=12)
ax21.legend(loc='upper right', frameon=False, fontsize=14)
plt.savefig(validation_plot_loc+'efield_first_three_halfwaves_'+plt_num+'.png', bbox_inches='tight', dpi=200)
#plt.savefig(validation_plot_loc+'efield_first_negative_'+plt_num+'.png', bbox_inches='tight', dpi=200)
plt.close()

# Spatially and temporally averaged values
full_t = time_range(time*1e6, full_pulse[0])
# spatially averaged
em_xavg = np.trapz(np.exp(em[:, block1]) * 6.022e23, x=xgrid[block1])/(xgrid[block1][-1] - xgrid[block1][0])
em_xtavg = np.trapz(em_xavg[full_t], x=time[full_t])/(time[full_t][-1] - time[full_t][0])
print(em_xtavg)
