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
    #exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_with_O2_V_m'+v+'_kV_10um.e')
    exopy = ExodusRead('/home/shane/projects/zapdos/problems/henry_interface/argon_water_henry_salt_out_07.e')

    plt_num = '001'

    emliq_density = exopy.get_vals('em_aq_density', data_type='element', block=1)
    H_density = exopy.get_vals('H_aq_density', data_type='element', block=1)
    OHm_density = exopy.get_vals('OHm_aq_density', data_type='element', block=1)
    Hp_density = exopy.get_vals('Hp_aq_density', data_type='element', block=1)
    OH_density = exopy.get_vals('OH__aqdensity', data_type='element', block=1)
    H2O2_density = exopy.get_vals('H2O2_aq_density', data_type='element', block=1)
    HO2_density = exopy.get_vals('HO2_aq_density', data_type='element', block=1)
    O2_density = exopy.get_vals('O2_aq_density', data_type='element', block=1)
    HO2m_density = exopy.get_vals('HO2m_aq_density', data_type='element', block=1)
    O3_density = exopy.get_vals('O3_aq_density', data_type='element', block=1)
    O2m_density = exopy.get_vals('O2m_aq_density', data_type='element', block=1)
    Om_density = exopy.get_vals('Om_aq_density', data_type='element', block=1)
    O3m_density = exopy.get_vals('O3m_aq_density', data_type='element', block=1)
    H2O_density = np.ones(np.shape(H_density)) * np.exp(10.92)*6.022e23

    k1 = 6e-8 / 6.022e23
    k2 = 3e7 / 6.022e23
    k3 = 1.1e7 / 6.022e23
    k4 = 3.5e3 / 6.022e23
    k5 = 1e-2 / 6.022e23
    k6 = 7e6 / 6.022e23
    k7 = 9e4 / 6.022e23
    k8 = 6e3 / 6.022e23
    k9 = 1.3e1 / 6.022e23
    k10 = 5.6e6 / 6.022e23
    k11 = 2.7e7 / 6.022e23
    k12 = 4.2e4 / 6.022e23
    k13 = 1.3e7 / 6.022e23
    k14 = 6e6 / 6.022e23
    k15 = 8e6 / 6.022e23
    k16 = 1.8e3 / 6.022e23
    k17 = 2.7e4 / 6.022e23
    k18 = 7.5e6 / 6.022e23
    k19 = 6e7 / 6.022e23
    k20 = 9e7 / 6.022e23
    k21 = 2.6e6 / 6.022e23
    k22 = 1.1e5 / 6.022e23
    k23 = 5e-4 / 6.022e23
    k24 = 5e-4 / 6.022e23
    k25 = 1.3e-4 / 6.022e23
    k26 = 9e7 / 6.022e23
    k27 = 1e2 / 6.022e23
    k28 = 1.6e5 / 6.022e23
    k29 = 5.3e6 / 6.022e23
    k30 = 3e9 / 6.022e23
        
    ts = -1
    rxn1 = k1 * emliq_density[ts, :] * H2Op_density[ts, :]
    rxn2 = -k2 * emliq_density[ts, :] * OH_density[ts, :]
    rxn3 = k3 * emliq_density[ts, :] * H2O2_density[ts, :]
    rxn4 = k4 * emliq_density[ts, :] * H2Om_density[ts, :]
    rxn5 = k5 * H_density[ts, :] * H2O_density[ts, :]
    rxn6 = -k6 * H_density[ts, :] * OH_density[ts, :]
    rxn7 = k7 * H_density[ts, :] * H2O2_density[ts, :]
    rxn8 = k8 * H2_density[ts, :] * _density[ts, :]
    rxn9 = k9 * O_density[ts, :] * _density[ts, :]
    rxn10 = -k10 * OH_density[ts, :] * OH_density[ts, :]
    rxn11 = -k11 * OH_density[ts, :] * Om_density[ts, :]
    rxn12 = -k12 * OH_density[ts, :] * H2_density[ts, :]
    rxn13 = -k13 * OH_density[ts, :] * OHm_density[ts, :]
    rxn14 = -k14 * OH_density[ts, :] * HO2_density[ts, :]
    rxn15 = -k15 * OH_density[ts, :] * O2m_density[ts, :]
    rxn16 = k16 * Om_density[ts, :] * H2O_density[ts, :]
    rxn17 = -k17 * OH_density[ts, :] * H2O2_density[ts, :]
    rxn18 = -k18 * OH_density[ts, :] * HO2m_density[ts, :]
    rxn19 = k19 * H3Op_density[ts, :] * OHm_density[ts, :]
    rxn20 = k20 * H_density[ts, :] * HO2_density[ts, :]
    rxn21 = -k21 * OH_density[ts, :] * O3m_density[ts, :]
    rxn22 = -k22 * OH_density[ts, :] * O3_density[ts, :]
    rxn23 = k23 * HO2_density[ts, :] * H2O2_density[ts, :]
    rxn24 = k24 * HO2_density[ts, :] * HO2m_density[ts, :]
    rxn25 = k25 * O2m_density[ts, :] * H2O2_density[ts, :]
    rxn26 = k26 * O3m_density[ts, :] * H_density[ts, :]
    rxn27 = k27 * HO3_density[ts, :] * O2_density[ts, :]
    rxn28 = k28 * O_density[ts, :] * H2O2_density[ts, :]
    rxn29 = k29 * O_density[ts, :] * HO2m_density[ts, :]
    rxn30 = k30 * O3_density[ts, :] * H2O2_density[ts, :]
                 #em_aq + H2Op_aq -> H_aq + OH_aq            : 6e-8 
                 #em_aq + OH_aq -> OHm_aq                    : 3e7
                 #em_aq + H2O2_aq -> OH_aq + OHm_aq                : 1.1e7 
                 #em_aq + HO2m_aq + H2O_aq -> OH_aq + OHm_aq + OHm_aq     : 3.5e3
                 #H_aq + H2O_aq -> H2_aq + OH_aq                   : 1e-2
                 #H_aq + OH_aq -> H2O_aq                        : 7e6
                 #H_aq + H2O2_aq -> OH_aq + H2O_aq                 : 9e4
                 #H2_aq + H2O2_aq -> H_aq + OH_aq + H2O_aq            : 6e3
                 #O_aq + H2O_aq -> OH_aq + OH_aq                   : 1.3e1
                 #OH_aq + OH_aq -> H2O2_aq    : 5.5e6
                 #OH_aq + Om_aq -> HO2m_aq    : 2e7
                 #OH_aq + H2_aq -> H_aq + H2O_aq     : 4.2e4
                 #OH_aq + OHm_aq -> Om_aq + H2O_aq   : 1.3e7
                 #OH_aq + HO2_aq -> H2O_aq + O2_aq   : 6e6
                 #OH_aq + O2m_aq -> OHm_aq + O2_aq     : 8e6
                 #Om_aq + H2O_aq -> OHm_aq + OH_aq     : 1.8e3
                 #OH_aq + H2O2_aq -> H2O_aq + HO2_aq     : 2.7e4
                 #OH_aq + HO2m_aq -> OHm_aq + HO2_aq     : 7.5e6
                 #H3Op_aq + OHm_aq -> H_aq + OH_aq + H2O_aq     : 6e7
                 #H_aq + HO2m_aq -> OHm_aq + OH_aq   : 9e7 
                 #OH_aq + O3m_aq ->  O3_aq + OHm_aq    : 2.6e6
                 #OH_aq + O3_aq -> HO2_aq + O2_aq          : 1.1e5
                 #HO2_aq + H2O2_aq -> OH_aq + O2_aq + H2O_aq    : 5e-4
                 #HO2_aq + HO2m_aq -> OHm_aq + OH_aq + O2_aq    : 5e-4
                 #O2m_aq + H2O2_aq -> OH_aq + O2_aq + OHm_aq            : 1.3e-4
                 #O3m_aq + H_aq -> O2_aq + OH_aq       : 9e7
                 #HO3_aq -> O2_aq + OH_aq           : 1e2
                 #O_aq + H2O2_aq -> OH_aq + HO2_aq     : 1.6e5
                 #O_aq + HO2m_aq -> OH_aq + O2m_aq     : 5.3e6
                 #O3_aq + H2O2_aq -> OH_aq + HO2_aq + O2_aq   : 3e9'

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
