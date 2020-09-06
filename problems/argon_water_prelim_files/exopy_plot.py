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
#voltage = ['1', '2', '3', '4', '5']
voltage = ['1']
num = 0
je_interface = np.empty(shape=(len(voltage)))
l_ne_aq = np.empty(shape=(len(voltage)), dtype='int')


f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,10))
for v in voltage:
    exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_model_prelim_V_m'+v+'_kV_10um.e')
    exopy2 = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_with_O2_V_m'+v+'_kV_10um.e')

    plt_num = '001'


    cmap = plt.get_cmap('jet')

    xgrid = exopy.exodus.variables['coordx'][:]

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

    block0 = exopy.split_block('x', 0)
    block1 = exopy.split_block('x', 1)

    cblock0 = block0[0:-1]
    cblock1 = block1[0:-1]

    # Get element data from nodes
    xcell = np.zeros(shape=(num_cells))
    for i in range(num_cells):
        xcell[i] = (xgrid[total_connect[i][1]-1] + xgrid[total_connect[i][0]-1])/2.0


    Eb0 = exopy.get_vals('Efield', data_type='element', block=0)
    Eb1 = exopy.get_vals('Efield', data_type='element', block=1)
    Efield = np.hstack((Eb0, Eb1))

    time = exopy.time[:]

    xb0, tb0 = np.meshgrid(xgrid[block0], exopy.time[:]*time_scale)
    xb1, tb1 = np.meshgrid(xgrid[block1], exopy.time[:]*time_scale)

    xc0, tb0 = np.meshgrid(xcell[block0[0:-2]], exopy.time[:]*time_scale)
    xc1, tc1 = np.meshgrid(xcell[block1[0:-2]], exopy.time[:]*time_scale)
    xc, tb = np.meshgrid(xcell, exopy.time[:]*time_scale)


    time_index = np.arange(len(exopy.time[:]))

    phi = exopy.get_vals('potential', data_type='node')[:]

    mean_en = exopy.get_vals('mean_en', data_type='node')[:]
    em = exopy.get_vals('em', data_type='node')[:]
    em_density = exopy.get_vals('em_density', data_type='element', block=0)
    emliq_density2 = exopy2.get_vals('emliq_density', data_type='element', block=1)
    emliq_density = exopy.get_vals('emliq_density', data_type='element', block=1)
    ne_aq = exopy.get_vals('emliq_density', data_type='element', block=1)
    Hp = exopy.get_vals('H+_density', data_type='element', block=1)
    OHm = exopy.get_vals('OHm_density', data_type='element', block=1)
    OH = exopy.get_vals('OH_density', data_type='element', block=1)
    je = exopy.get_vals('Current_em', data_type='element', block=0)
    je_aq = exopy.get_vals('Current_emliq', data_type='element', block=1)
    johm = exopy.get_vals('Current_OHm', data_type='element', block=1)
    #je_interface[num] = je[-1,-1]
    je_interface[num] = je_aq[-1,0]

    #l_ne_aq[num] = np.where(ne_aq[-1,:] <= 0.1*ne_aq[-1,0])[:][0][0]
    #l_ne_aq[num] = np.where(ne_aq[-1, :] <= 1e20)[:][0][0]
    l_ne_aq[num] = np.where(je_aq[-1, :] >= -0.1)[:][0][0]

    plt.semilogy(xcell[cblock1], emliq_density[-1,:])
    plt.semilogy(xcell[cblock1], emliq_density2[-1, :])

    # color cycle:
    color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    #plt.plot(xcell[cblock1], OHm[-1,:], label=v+' kV')
    #plt.semilogy(xcell[cblock1], OH[-1,:], label=v+' kV')
    #plt.plot(xcell[cblock0], Eb0[-1,:], color_cycle[num])
    #plt.plot(xcell[cblock1], Eb1[-1,:], color_cycle[num], label=v+' kV')
    #plt.plot(xcell[cblock0], je[-1,:], label=v+' kV')
    #plt.plot(xcell[cblock0], je[-1,:], color_cycle[num])
    #plt.plot(xcell[cblock1], je_aq[-1, :], color_cycle[num])
    #plt.plot(xcell[cblock1], johm[-1, :], color_cycle[num], linestyle='--')
    #if (v == '1old'):
    #    #plt.semilogy(xcell[cblock1], ne_aq[-1,:], color_cycle[num-1], linestyle='--', label=v+' kV')
    #    plt.semilogy(xcell[cblock1], ne_aq[-1,:], color_cycle[num-1], label=v+' kV')
    #else:
    #    plt.semilogy(xcell[cblock1], ne_aq[-1,:], color_cycle[num], label=v+' kV')
    #plt.semilogy(xcell[cblock1], ne_aq[-1,:]*ne_aq[-1,:]*5.5e6/6.022e23)

    #for x in range(len(liquid_species)):
    #    var = exopy.get_vals(liquid_species[x]+'_density', data_type='element', block=1)
    #    plt.semilogy(xcell[cblock1], var[-1, :], label=liquid_species[x])
    #plt.legend()
    #plt.show()
    #exit()

    #plt.plot(xgrid, phi[-1,:], label=v+' kV')
    #two_region_plot(xgrid[block0], xgrid[block1], phi[-1,block0], phi[-1,block1])
    #two_region_plot(xcell[cblock0], xcell[cblock1], Eb0[-1,:], Eb1[-1,:])



    #ax1.plot(xcell[cblock0]*1e3, em_density[-1, :], 'k')
    #ax1.set_xlim([0, 1])
    ##ax1.set_xlim([0.8, 1])
    ##ax1.set_ylim([-0.2e7, 0.1e7])

    #ax2.plot(xcell[cblock1]*1e3, emliq_density[-1, :], 'k')
    #ax2.set_xlim([1, 1.01])
    ##ax2.set_xlim([1, 1.0001])

    #xticks = ax1.xaxis.get_major_ticks()
    #xticks[-1].label1.set_visible(False)

    ##ax2.axes.tick_params(axis='y', which='both', bottom=False, top=False)
    #ax2.yaxis.set_visible(False)

    #ax1.xaxis.set_major_formatter(mtick.FormatStrFormatter('%1.2f'))
    #ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%1.2f'))
    ##xfmt = mtick.ScalarFormatterForceFormat()
    ##xfmt.set_powerlimits((0,0))
    ##ax1.xaxis.set_major_formatter(xfmt)

    #f.text(0.5, 0.04, 'X [mm]', ha='center')
    #f.subplots_adjust(wspace=0)

    num += 1

#plt.plot(abs(je_interface), xcell[l_ne_aq], '-s')
plt.legend()
plt.show()


#plt.xlabel('Current Density (A m$^2$)', fontsize=16)
#plt.ylabel('Penetration Depth (m)', fontsize=16)
#plt.savefig('penetration_depth_vs_current_density.png', dpi=200, bbox_inches='tight')
#plt.close()
