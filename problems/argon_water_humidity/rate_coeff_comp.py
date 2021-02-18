import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

file01 = 'bolsig_files/'
file05 = 'bolsig_files_05/'
file10 = 'bolsig_files_10/'

#reaction = 'C1_Ar_Elastic'
reaction = 'C60_H2O_Elastic'
#reaction = 'C67_H2O_Ionization_13.50_eV'

data01 = pd.read_csv(file01+reaction, sep=' ', names=['x', 'y'])
data05 = pd.read_csv(file05+reaction, sep=' ', names=['x', 'y'])
data10 = pd.read_csv(file10+reaction, sep=' ', names=['x', 'y'])

plt.plot(data01['x'], data01['y'], 'o-')
plt.plot(data05['x'], data05['y'], 'o-')
plt.plot(data10['x'], data10['y'], 'o-')
plt.show()
