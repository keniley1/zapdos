from exodus_read_base2 import ExodusRead
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick

voltage = ['1', '2', '3', '4', '5']
#voltage = ['1', '2', '3', '4']
num = 0
je_interface = np.empty(shape=(len(voltage)))
l_ne_aq = np.empty(shape=(len(voltage)), dtype='int')


gas_species = ['Arp', 'em']
liquid_species = ['emliq', 'OHm', 'Om', 'O2m', 'O3m', 'HO2m', 'H+', 
                  'O', 'O2_1', 'O3', 'H', 'H2', 'HO2', 'OH', 'H2O2']

plot_species = 'O2'

for v in voltage:
    #exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_model_prelim_V_m'+v+'_kV_10um.e')
    exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_with_O2_V_m'+v+'_kV_10um.e')

    plt_num = '001'


    cmap = plt.get_cmap('jet')

    xgrid = exopy.exodus.variables['coordx'][:]

    time_scale = 1

    dt_start = 0
    dt_end = -1


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


    xc1, tc1 = np.meshgrid(xcell[block1[0:-2]], exopy.time[:]*time_scale)


    density = exopy.get_vals(plot_species+'_density', data_type='element', block=1)


    #plt.semilogy(xcell[cblock1]*1e3, emliq_density[-1, :], label=v+' kV')
    plt.semilogy(xcell[cblock1]*1e3, density[-1, :], label=v+' kV')






plt.xlabel('X (mm)', fontsize=14)
plt.ylabel('Density (m$^{-3}$)', fontsize=14)
#plt.show()
plt.legend(loc='best', frameon=False)

#plt.title('H$_2$O$_2$', fontsize=16)
#plt.savefig('plots/H2O2_varying_voltage.png', dpi=200, bbox_inches='tight')
#plt.close()

#plt.title('H$_2$', fontsize=16)
#plt.savefig('plots/H2_varying_voltage.png', dpi=200, bbox_inches='tight')
#plt.close()


#plt.title('OH$^-$', fontsize=16)
#plt.savefig('plots/OHm_varying_voltage.png', dpi=200, bbox_inches='tight')
#plt.close()

#plt.title('e$_{aq}$', fontsize=16)
#plt.savefig('plots/solvated_electrons_varying_voltage.png', dpi=200, bbox_inches='tight')
#plt.close()

plt.title('O$_2$', fontsize=16)
plt.savefig('plots/O2_varying_voltage.png', dpi=200, bbox_inches='tight')
plt.close()
