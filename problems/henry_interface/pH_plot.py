from exodus_read_base2 import ExodusRead
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick


voltage = ['1', '2', '3', '4', '5']
#voltage = ['1']
num = 0
pH_v = np.empty(shape=(len(voltage)))
pOH_v = np.empty(shape=(len(voltage)))
je_interface = np.empty(shape=(len(voltage)))

f, ax2 = plt.subplots(1, 1, figsize=(8,8))
for v in voltage:
    exopy = ExodusRead('/home/shane/projects/zapdos/problems/argon_water_prelim_files/argon_water_gopalakrishnan_m'+v+'_kV_1um.e')

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

    x0 = xgrid[block0]
    x1 = xgrid[block1]

    time = exopy.time[:]

    xb0, tb0 = np.meshgrid(xgrid[block0], exopy.time[:]*time_scale)
    xb1, tb1 = np.meshgrid(xgrid[block1], exopy.time[:]*time_scale)

    time_index = np.arange(len(exopy.time[:]))

    Hp_density = np.exp(exopy.get_vals('H+', data_type='node')[:,block1])*6.022e23
    H3Op_density = np.exp(exopy.get_vals('H3O+', data_type='node')[:,block1])*6.022e23
    OHm_density = np.exp(exopy.get_vals('OH-', data_type='node')[:,block1])*6.022e23
    Hp_molarity = (Hp_density + H3Op_density)/6.022e23 * 1e-3
    OHm_molarity = OHm_density/6.022e23 * 1e-3
    pH = -np.log10(Hp_molarity)
    pOH = -np.log10(OHm_molarity)

    emliq_density = exopy.get_vals('emliq', data_type='node')[:,block1]
    emliq_density = np.exp(emliq_density)*6.022e23
    potential = exopy.get_vals('potential', data_type='node')[:,block1]

    mu = 0.000173913 / 1000.
    diff = 4.5e-9 
    je_aq = -(-mu * -np.gradient(potential[-1,:], x1) * emliq_density[-1,:] - diff * np.gradient(emliq_density[-1,:], x1)) * 1.602e-19
    je_interface[num] = je_aq[0]

    # color cycle:
    color_cycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    #ax2.semilogy(x1*1e3, Hp_density[-1, :], label='H$^+$')
    #ax2.semilogy(x1*1e3, H3Op_density[-1, :], label='H3O$^+$')
    #plt.plot(x1*1e3, pH[-1, :], label='H$^+$')
    #plt.semilogx(time, pH[:, -1], label='V = -'+str(v)+' kV')
    #ax2.set_xlim([1, 1.001])
    #plt.plot(v, pH[-1, 0], 'ko', label='V = -'+str(v)+' kV')
    #plt.plot(v, pOH[-1, 0], 'ro')
    #ax2.legend(loc='best', frameon=False, ncol=3)

    #ax2.axes.tick_params(axis='y', which='both', bottom=False, top=False)

    #ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%1.3f'))
    #ax2.tick_params(axis='both', labelsize=15)

    #f.subplots_adjust(wspace=0)
    pH_v[num] = pH[-1,0]
    pOH_v[num] = pOH[-1,0]

    num += 1

#plt.plot(voltage, pH_v, 'ko-', label='pH')
plt.plot(abs(je_interface), pH_v, 'ko-', label='pH')
plt.plot()
ax2.tick_params(axis='both', labelsize=15)
plt.xlabel('Electron Current Density, j$_e$ (A m$^{-2}$)', fontsize=18)
plt.ylabel('pH', fontsize=18)
plt.title('pH at Plasma-Water Interface', fontsize=20)
#plt.show()
#exit()
plt.savefig('plots/pH_interface.png', dpi=200, bbox_inches='tight')
plt.close()
exit()
