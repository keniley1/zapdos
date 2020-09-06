from exodus_read_base2 import ExodusRead
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick

class ScalarFormatterForceFormat(mtick.ScalarFormatter):
    def _set_format(self,vmin,vmax):
        self.format = "%1.1f"

def two_region_plot(x1, x2, y1, y2):
    #fig = plt.figure(figsize = (10,10))
    #gs1 = gridspec.GridSpec(1,2)
    #gs1.update(wspace=0.025, hspace=0.025)


    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,10))
    ax1.plot(x1*1e3, y1)
    ax1.set_xlim([0, 1])
    #ax1.set_xlim([0.8, 1])
    #ax1.set_ylim([-0.2e7, 0.1e7])

    ax2.plot(x2*1e3, y2)
    ax2.set_xlim([1, 1.01])
    #ax2.set_xlim([1, 1.0001])

    xticks = ax1.xaxis.get_major_ticks()
    xticks[-1].label1.set_visible(False)

    #ax2.axes.tick_params(axis='y', which='both', bottom=False, top=False)
    ax2.yaxis.set_visible(False)

    ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%1.2f'))
    ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%1.2f'))
    #xfmt = mtick.ScalarFormatterForceFormat()
    #xfmt.set_powerlimits((0,0))
    #ax1.xaxis.set_major_formatter(xfmt)

    f.text(0.5, 0.04, 'X [mm]', ha='center')
    f.subplots_adjust(wspace=0)
    #plt.show()
    

def extruded_plot(xmesh, ymesh, variable, xlabel, ylabel, title='', dpi=300, save=True, save_name='none', show=False, time_data=0):
    if (save and show):
        print("ERROR: cannot save and show plot! Choose one.")
        exit()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    m0 = ax.pcolormesh(xmesh, ymesh, variable, cmap='jet', norm=colors.LogNorm(vmin=variable.min(), vmax=variable.max()))
    #m0.set_clim(-3.6e6, 3.6e6)
    #m0.set_clim(0, 1e18)
    fig.colorbar(m0)
    plt.xlabel('Gap distance (m)')
    plt.ylabel('Time ($\mu$s)')
    plt.title(title)
    plt.savefig(save_name+'.png', dpi=dpi, bbox_inches='tight')
    plt.close()

def line_plot_peaks(x, y, time, include_potential=False, potential=0): 
    if (include_potential and not potential):
        print("ERROR: If include_potential == True, the time array and potential array needs to be included as well!")
        exit()

    fig = plt.figure(figsize=(8,6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1,3])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    min_indices = np.where((np.r_[True, potential[1:] < potential[:-1]] & np.r_[potential[:-1] < potential[1:], True])==True)[0]
    max_indices = np.where((np.r_[True, potential[1:] > potential[:-1]] & np.r_[potential[:-1] > potential[1:], True])==True)[0]
    vv = np.linspace(-10, 10, 10)
    min1 = np.zeros(shape=(10,)) + time[min_indices[0]] 
    min2 = np.zeros(shape=(10,)) + time[min_indices[1]]
    max1 = np.zeros(shape=(10,)) + time[max_indices[0]]
    max2 = np.zeros(shape=(10,)) + time[max_indices[1]]
    ax0.plot(time, Vw, 'k')
    ax0.plot(min1, vv, c='#1f77b4')
    ax0.plot(min2, vv, c='#2ca02c')
    ax0.plot(max1, vv, c='#ff7f0e')
    ax0.plot(max2, vv, c='#d62728')
    ax0.set_xlim([0, 8e-6])
    ax0.set_ylim([-10, 10])
    ax0.set_ylabel('Potential [kV]', fontsize=14)
    for i in range(len(min_indices)):
        ax1.semilogy(x, y[:,min_indices])
        ax1.semilogy(x, y[:,max_indices])
    plt.show()


def line_plot(x, y, time, include_potential=False, potential=0): 
    print("Not there yet...")
    exit()


#voltage = ['1', '1old', '2', '3', '4', '5']
voltage = ['1']
num = 0
je_interface = np.empty(shape=(len(voltage)))
l_ne_aq = np.empty(shape=(len(voltage)), dtype='int')


f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(12,10))
for v in voltage:
    exopy = ExodusRead('/home/shane/projects/zapdos/problems/henry_interface/argon_water_henry_out_01.e')

    plt_num = '001'


    cmap = plt.get_cmap('jet')

    xgrid = exopy.exodus.variables['coordx'][:]

    time_scale = 1

    dt_start = 0
    dt_end = -1

    gas_species = ['OH']
    liquid_species = ['OH_aq']

    #print(exopy.exodus.variables.keys())

    cell_block0 = exopy.exodus.variables['connect1'][:]
    cell_block1 = exopy.exodus.variables['connect2'][:]

    total_connect = np.vstack((cell_block0, cell_block1))
    num_cells = len(total_connect)

    block0 = exopy.split_block('x', 0)
    block1 = exopy.split_block('x', 1)

    x0 = xgrid[block0]
    x1 = xgrid[block1]

    time = exopy.time[:]

    xb0, tb0 = np.meshgrid(xgrid[block0], exopy.time[:]*time_scale)
    xb1, tb1 = np.meshgrid(xgrid[block1], exopy.time[:]*time_scale)

    time_index = np.arange(len(exopy.time[:]))

    OH_density = np.exp(exopy.get_vals('OH', data_type='node')[:,block0])*6.022e23

    emliq_density = np.exp(exopy.get_vals('emliq', data_type='node')[:,block1])*6.022e23
    OHm_density = np.exp(exopy.get_vals('OH-', data_type='node')[:,block1])*6.022e23
    OH_aq_density = np.exp(exopy.get_vals('OH_aq', data_type='node')[:,block1])*6.022e23
    #jemliq = exopy.get_vals('Current_emliq', data_type='element', block=1)
    #jem = exopy.get_vals('Current_em', data_type='element', block=0)

    # color cycle:
    color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    ax1.semilogy(x0*1e3, OH_density[-1, :], 'r', linewidth=3, label='OH$_g$')
    ax1.set_ylabel('Density (m$^{-3}$)', fontsize=16)
    ax1.set_xlim([0, 1])
    #ax1.set_xlim([0.8, 1])
    ax1.set_ylim([1e14, 1e17])
    ax1.legend(loc='best', frameon=False)

    #ax2.semilogy(x1*1e3, emliq_density[-1, :], 'k', label='e$_{aq}$')
    ax2.semilogy(x1*1e3, OH_aq_density[-1, :], 'r', linewidth=3, label='OH$_{aq}$')
    ax2.set_xlim([1, 1.001])
    ax2.legend(loc='best', frameon=False, ncol=3)
    #ax2.set_xlim([1, 1.0001])

    xticks = ax1.xaxis.get_major_ticks()
    xticks[-1].label1.set_visible(False)

    #ax2.axes.tick_params(axis='y', which='both', bottom=False, top=False)
    ax2.yaxis.set_visible(False)

    ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%1.3f'))
    ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%1.3f'))
    #xfmt = mtick.ScalarFormatterForceFormat()
    #xfmt.set_powerlimits((0,0))
    #ax1.xaxis.set_major_formatter(xfmt)
    ax1.tick_params(axis='both', labelsize=15)
    ax2.tick_params(axis='both', labelsize=15)

    f.text(0.5, 0.04, 'X [mm]', ha='center', fontsize=18)
    f.text(0.5, 0.9, 'V$_{DC}$ = -1 kV', ha='center', fontsize=18)
    f.subplots_adjust(wspace=0)

    num += 1

plt.plot()
#plt.show()
#exit()
plt.savefig('plots/oh_interface.png', dpi=200, bbox_inches='tight')
plt.close()
exit()

# Interface layer time!
plt.semilogy(xcell[cblock1]*1e3, emliq_density[-1, :], 'k', label='e$_{aq}$')
plt.semilogy(xcell[cblock1]*1e3, OHm_density[-1, :], label='OH$^-$')
plt.semilogy(xcell[cblock1]*1e3, Om_density[-1, :], label='O$^-$')
plt.semilogy(xcell[cblock1]*1e3, O2m_density[-1, :], label='O$_2^-$')
plt.semilogy(xcell[cblock1]*1e3, O3m_density[-1, :], label='O$_3^-$')
plt.semilogy(xcell[cblock1]*1e3, HO2m_density[-1, :], label='HO$_2^-$')
plt.semilogy(xcell[cblock1]*1e3, Hp_density[-1, :], label='H$^+$')
plt.semilogy(xcell[cblock1]*1e3, O_density[-1, :], label='O')
plt.semilogy(xcell[cblock1]*1e3, O2_1_density[-1, :], label='O$_2^1$')
plt.semilogy(xcell[cblock1]*1e3, O3_density[-1, :], label='O$_3$')
plt.semilogy(xcell[cblock1]*1e3, H_density[-1, :], label='H')
plt.semilogy(xcell[cblock1]*1e3, H2_density[-1, :], label='H$_2$')
plt.semilogy(xcell[cblock1]*1e3, HO2_density[-1, :], label='HO$_2$')
plt.semilogy(xcell[cblock1]*1e3, OH_density[-1, :], label='OH')
plt.semilogy(xcell[cblock1]*1e3, H2O2_density[-1, :], label='H$_2$O$_2$')
plt.axis([1, 1.001, 1e14, 1e23])
plt.xlabel('X [mm]', fontsize=12)
plt.ylabel('Density (m$^{-3}$)', fontsize=12)
plt.savefig('plots/interface_zoom.png', dpi=200, bbox_inches='tight')
plt.close()


#plt.xlabel('Current Density (A m$^2$)', fontsize=16)
#plt.ylabel('Penetration Depth (m)', fontsize=16)
#plt.savefig('penetration_depth_vs_current_density.png', dpi=200, bbox_inches='tight')
#plt.close()
