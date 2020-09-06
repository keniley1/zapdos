from exodus_read_base2 import ExodusRead
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit

voltage = ['1', '2', '3', '4', '5']
#voltage = ['1', '2', '3', '4']
num = 0
je_interface = np.empty(shape=(len(voltage)))
l_rumbach = np.empty(shape=(len(voltage)))
l_simulation = np.empty(shape=(len(voltage)))
l_test = np.empty(shape=(len(voltage)))
l_simulation2 = np.empty(shape=(len(voltage)))
n0_rumbach = np.empty(shape=(len(voltage)))
n0_simulation = np.empty(shape=(len(voltage)))

def func(x, a, b):
    return(a * np.exp(-x/b))


for v in voltage:
    #exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_model_prelim_V_m'+v+'_kV_10um.e')
    exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_with_O2_V_m'+v+'_kV_10um.e')

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

    xcell_liq = xcell[cblock1]

    xb0, tb0 = np.meshgrid(xgrid[block0], exopy.time[:]*time_scale)
    xb1, tb1 = np.meshgrid(xgrid[block1], exopy.time[:]*time_scale)

    xc0, tb0 = np.meshgrid(xcell[block0[0:-2]], exopy.time[:]*time_scale)
    xc1, tc1 = np.meshgrid(xcell[block1[0:-2]], exopy.time[:]*time_scale)
    xc, tb = np.meshgrid(xcell, exopy.time[:]*time_scale)

    emliq_density = exopy.get_vals('emliq_density', data_type='element', block=1)
    em_density = exopy.get_vals('em_density', data_type='element', block=0)

    popt, pcov = curve_fit(func, xcell_liq-1e-3, emliq_density[-1,:], p0=[emliq_density[-1,0], 1e-8], method='trf')

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

    l_simulation2[num] = xcell_liq[-1]
    for xx in range(len(johm[-1,:])):
        if (abs(johm[-1,xx]) >= abs(je_aq[-1,xx])):
            l_simulation2[num] = xcell_liq[xx] - 1e-3
            break

    # Faraday constant 
    F = 9.65e4
    
    # Interfacial concentration
    n0_rumbach[num] = (abs(je[-1,-1])**2.0 / (2.0*5.5e9*4.5e-9*(F**2.0)))**(1.0/3.0)
    n0_simulation[num] = emliq_density[-1,0]/6.022e23

    # Penetration depth
    l_rumbach[num] = (F*4.5e-9/(2.0*5.5e9*-1*je[-1,-1]))**(1.0/3.0)

    l_test[num] = popt[1]

    # Exponential fit
    #n_test = n0_simulation[num] * 6.022e23 * np.exp(-(xcell_liq-1e-3) / l_rumbach[num])
    n_test = func(xcell_liq-1e-3, *popt)
    plt.semilogy(xcell_liq, emliq_density[-1,:])
    plt.semilogy(xcell_liq, n_test)
    #plt.xlim([1e-3, 1.001e-3])
    #plt.ylim([1e20, 1e23])
    plt.show()

    num += 1

plt.plot(abs(je_interface), l_rumbach*1e6, '-s', label='Rumbach Estimate')
plt.plot(abs(je_interface), l_simulation*1e6, '-s', label='Zapdos-Crane')
plt.plot(abs(je_interface), l_test*1e6, '-s', label='Test')
plt.ylabel('Penetration Depth, l ($\mu$m)', fontsize=14)
plt.xlabel('Electron Current Density, j$_e$ (A m$^{-2}$)', fontsize=14)
plt.legend(frameon=False)
#plt.savefig('plots/rumbach_penetration_depth_comparison.png', dpi=200, bbox_inches='tight')
#plt.close()
plt.show()
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
