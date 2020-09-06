from exodus_read_base2 import ExodusRead
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit

voltage = ['1', '2', '3', '4', '5']
num = 0
je_interface = np.empty(shape=(len(voltage)))
l_rumbach = np.empty(shape=(len(voltage)))
l_simulation = np.empty(shape=(len(voltage)))
n0_rumbach = np.empty(shape=(len(voltage)))
n0_simulation = np.empty(shape=(len(voltage)))

#reactions: emliq + H2O -> H + OHm     : 1.9e-2 
#                 H + OHm -> H2O + emliq     : 2.2e4
#                 H -> emliq + H+            : 3.9
#                 emliq + H+ -> H            : 2.3e7
#                 emliq + OH -> OHm          : 3.0e7
#                 emliq + H2O2 -> OH + OHm   : 1.1e7
#                 emliq + HO2 -> HO2m        : 2.0e7
#                 emliq + O2 -> O2m          : 1.9e7
#                 emliq + HO2m -> Om + OHm   : 3.5e7
#                 emliq + O3 -> O3m          : 3.6e7
#                 emliq + O2m -> HO2m + OHm          : 1.3e7
#                 emliq + H -> H2 + OHm              : 2.5e7
#                 emliq + Om -> OHm + OHm            : 2.2e7
#                 emliq + O3m -> O2 + OHm + OHm      : 1.6e7
#                 emliq + emliq -> H2 + OHm + OHm    : 5.5e6

def func(x, a, b):
    return(a * np.exp(-x/b))


for v in voltage:
    #exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_model_prelim_V_m'+v+'_kV_10um.e')
    exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_with_O2_V_m'+v+'_kV_10um.e')

    plt_num = '001'

    emliq_density = exopy.get_vals('emliq_density', data_type='element', block=1)
    em_density = exopy.get_vals('em_density', data_type='element', block=0)
    H_density = exopy.get_vals('H_density', data_type='element', block=1)
    OHm_density = exopy.get_vals('OHm_density', data_type='element', block=1)
    Hp_density = exopy.get_vals('H+_density', data_type='element', block=1)
    OH_density = exopy.get_vals('OH_density', data_type='element', block=1)
    H2O2_density = exopy.get_vals('H2O2_density', data_type='element', block=1)
    HO2_density = exopy.get_vals('HO2_density', data_type='element', block=1)
    O2_density = exopy.get_vals('O2_density', data_type='element', block=1)
    HO2m_density = exopy.get_vals('HO2m_density', data_type='element', block=1)
    O3_density = exopy.get_vals('O3_density', data_type='element', block=1)
    O2m_density = exopy.get_vals('O2m_density', data_type='element', block=1)
    Om_density = exopy.get_vals('Om_density', data_type='element', block=1)
    O3m_density = exopy.get_vals('O3m_density', data_type='element', block=1)
    H2O_density = np.ones(np.shape(H_density)) * np.exp(10.92)*6.022e23

    k1 = 1.9e-2 / 6.022e23
    k2 = 2.2e4 / 6.022e23
    k3 = 3.9 / 6.022e23
    k4 = 2.3e7 / 6.022e23
    k5 = 3e7 / 6.022e23
    k6 = 1.1e7 / 6.022e23
    k7 = 2e7 / 6.022e23
    k8 = 1.9e7 / 6.022e23
    k9 = 3.5e7 / 6.022e23
    k10 = 3.6e7 / 6.022e23
    k11 = 1.3e7 / 6.022e23
    k12 = 2.5e7 / 6.022e23
    k13 = 2.2e7 / 6.022e23
    k14 = 1.6e7 / 6.022e23
    k15 = 5.5e6 / 6.022e23
        
    ts = -1
    rxn1 = k1 * emliq_density[ts, :] * H2O_density[ts, :]
    rxn2 = k2 * H_density[ts, :] * OHm_density[ts, :]
    rxn3 = k3 * H_density[ts, :] 
    rxn4 = k4 * emliq_density[ts, :] * Hp_density[ts, :]
    rxn5 = k5 * emliq_density[ts, :] * OH_density[ts, :]
    rxn6 = k6 * emliq_density[ts, :] * H2O2_density[ts, :]
    rxn7 = k7 * emliq_density[ts, :] * HO2_density[ts, :]
    rxn8 = k8 * emliq_density[ts, :] * O2_density[ts, :]
    rxn9 = k9 * emliq_density[ts, :] * HO2m_density[ts, :]
    rxn10 = k10 * emliq_density[ts, :] * O3_density[ts, :]
    rxn11 = k11 * emliq_density[ts, :] * O2m_density[ts, :]
    rxn12 = k12 * emliq_density[ts, :] * H_density[ts, :]
    rxn13 = k13 * emliq_density[ts, :] * Om_density[ts, :]
    rxn14 = k14 * emliq_density[ts, :] * O3m_density[ts, :]
    rxn15 = k15 * emliq_density[ts, :] * emliq_density[ts, :]

    cmap = plt.get_cmap('jet')

    xgrid = exopy.exodus.variables['coordx'][:]

    time_scale = 1

    dt_start = 0
    dt_end = -1

    gas_species = ['Arp', 'em']
    liquid_species = ['emliq', 'OHm', 'Om', 'O2m', 'O3m', 'HO2m', 'H+', 
                  'O', 'O2_1', 'O3', 'H', 'H2', 'HO2', 'OH', 'H2O2']

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

    xcell_liq = xcell[cblock1]

    color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

#reactions: emliq + H2O -> H + OHm     : 1.9e-2 
#                 H + OHm -> H2O + emliq     : 2.2e4
#                 H -> emliq + H+            : 3.9
#                 emliq + H+ -> H            : 2.3e7
#                 emliq + OH -> OHm          : 3.0e7
#                 emliq + H2O2 -> OH + OHm   : 1.1e7
#                 emliq + HO2 -> HO2m        : 2.0e7
#                 emliq + O2 -> O2m          : 1.9e7
#                 emliq + HO2m -> Om + OHm   : 3.5e7
#                 emliq + O3 -> O3m          : 3.6e7
#                 emliq + O2m -> HO2m + OHm          : 1.3e7
#                 emliq + H -> H2 + OHm              : 2.5e7
#                 emliq + Om -> OHm + OHm            : 2.2e7
#                 emliq + O3m -> O2 + OHm + OHm      : 1.6e7
#                 emliq + emliq -> H2 + OHm + OHm    : 5.5e6

    plt.semilogy(xcell_liq, rxn1,        color=color_cycle[0], label='H$_2$O')
    plt.semilogy(xcell_liq, rxn4,        color=color_cycle[1], label='H$^+$')
    plt.semilogy(xcell_liq, rxn5,        color=color_cycle[2], label='OH')
    plt.semilogy(xcell_liq, rxn6,        color=color_cycle[3], label='H$_2$O$_2$')
    plt.semilogy(xcell_liq, rxn7,        color=color_cycle[4], label='HO$_2$')
    plt.semilogy(xcell_liq, rxn8,        color=color_cycle[5], label='O$_2$')
    plt.semilogy(xcell_liq, rxn9,        color=color_cycle[6], label='HO$_2^-$')
    plt.semilogy(xcell_liq, rxn10,       color=color_cycle[7], label='O$_3$')
    plt.semilogy(xcell_liq, rxn11,       color=color_cycle[8], label='O$_2^-$')
    plt.semilogy(xcell_liq, rxn12,       color=color_cycle[9], label='H')
    plt.semilogy(xcell_liq, rxn13, '--', color=color_cycle[0], label='O$^-$')
    plt.semilogy(xcell_liq, rxn14, '--', color=color_cycle[1], label='O$_3^-$')
    plt.semilogy(xcell_liq, rxn15, '--', color=color_cycle[2], label='e$_{(aq)}$')
    # These two are the two sources. Probably irrelevant.
    plt.semilogy(xcell_liq, rxn2, '--', color=color_cycle[3], label='H + OH$^-$ (sink)')
    plt.semilogy(xcell_liq, rxn3, '--', color=color_cycle[4], label='H (sink)')
    plt.legend(loc='best')
    plt.show()
    exit()
    plt.close()

    # Rxn15 is clearly the highest by ~10x
    # Let's lump the rest together and see what it looks like...
    plt.semilogy(xcell_liq, rxn15, label='e$_{(aq)}$')
    plt.semilogy(xcell_liq, rxn1 + rxn2 + rxn3 + rxn4 + rxn5 + rxn6 + rxn7 + rxn8 + rxn9 + rxn10 + rxn11 + rxn12 + rxn13 + rxn14, label='all')
    plt.show()

    xb0, tb0 = np.meshgrid(xgrid[block0], exopy.time[:]*time_scale)
    xb1, tb1 = np.meshgrid(xgrid[block1], exopy.time[:]*time_scale)

    xc0, tb0 = np.meshgrid(xcell[block0[0:-2]], exopy.time[:]*time_scale)
    xc1, tc1 = np.meshgrid(xcell[block1[0:-2]], exopy.time[:]*time_scale)
    xc, tb = np.meshgrid(xcell, exopy.time[:]*time_scale)

    popt, pcov = curve_fit(func, xcell_liq, emliq_density[-1,:], bounds=(0, 0.1))

    time_index = np.arange(len(exopy.time[:]))

    phi = exopy.get_vals('potential', data_type='node')[:]

    je = exopy.get_vals('Current_em', data_type='element', block=0)
    je_aq = exopy.get_vals('Current_emliq', data_type='element', block=1)
    johm = exopy.get_vals('Current_OHm', data_type='element', block=1)
    #je_interface[num] = je[-1,-1]
    je_interface[num] = je_aq[-1,0]

    l_simulation[num] = xcell_liq[-1]
    #for xx in range(len(johm[-1,:])):
    try:
        temp = np.where(emliq_density[-1,:] <= emliq_density[-1,0]*0.03)[:][0][0]
    except:
        temp = -1
    if (xcell_liq[temp] < l_simulation[num]):
        l_simulation[num] = xcell_liq[temp] - 1e-3

    # Faraday constant 
    F = 9.65e4
    
    # Interfacial concentration
    n0_rumbach[num] = (abs(je[-1,-1])**2.0 / (2.0*5.5e9*4.5e-9*(F**2.0)))**(1.0/3.0)
    n0_simulation[num] = emliq_density[-1,0]/6.022e23

    # Penetration depth
    l_rumbach[num] = (F*4.5e-9/(2.0*5.5e9*-1*je[-1,-1]))**(1.0/3.0)

    num += 1


#plt.plot(abs(je_interface), l_rumbach*1e6, '-s', label='Rumbach Estimate')
#plt.plot(abs(je_interface), l_simulation*1e6, '-s', label='Zapdos-Crane')
#plt.ylabel('Penetration Depth, l ($\mu$m)', fontsize=14)
#plt.xlabel('Electron Current Density, j$_e$ (A m$^{-2}$)', fontsize=14)
#plt.legend(frameon=False)
#plt.savefig('plots/rumbach_penetration_depth_comparison.png', dpi=200, bbox_inches='tight')
#plt.close()
#plt.show()
exit()

plt.plot(abs(je_interface), n0_rumbach, '-s', label='Rumbach Estimate')
plt.plot(abs(je_interface), n0_simulation, '-s', label='Zapdos-Crane')
plt.ylabel('Interfacial Concentration, n$_0$ (mM)', fontsize=14)
plt.xlabel('Electron Current Density, j$_e$ (A m$^{-2}$)', fontsize=14)
plt.legend(frameon=False)
plt.savefig('plots/rumbach_interfacial_concentration_comparison.png', dpi=200, bbox_inches='tight')
plt.close()

exit()

#plt.xlabel('Current Density (A m$^2$)', fontsize=16)
#plt.ylabel('Penetration Depth (m)', fontsize=16)
#plt.savefig('penetration_depth_vs_current_density.png', dpi=200, bbox_inches='tight')
#plt.close()
