import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data0 = pd.read_csv('gas_temp_out0.csv')
data1 = pd.read_csv('gas_temp_out1.csv')

data3500 = pd.read_csv('gas_temp_3500.csv')

xg = data0['x_node']
Tg = data0['Tg']

Tw = data1['Tw']
xw = data1['x_node']

h2o = data0['H2O']
h2o_350 = data3500['H2O']

#fig, (ax1,ax2) = plt.subplots(2,1, figsize=(8,10))

#ax1 = plt.plot(xg, Tg)
#ax2 = plt.plot(xw, Tw)

# eqn
T = 300
p = 0.61094*exp(17.625*T/(T + 243.04))

#plt.plot(xg, np.exp(h2o)*6.022e23)
#plt.plot(xg, np.exp(h2o_350)*6.022e23)

plt.show()
