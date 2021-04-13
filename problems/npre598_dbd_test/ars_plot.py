import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

nt = 2000
time = np.zeros(shape=(nt))

output_loc = 'output_data/'
for i in range(nt):
    data = pd.read_csv(output_loc+'node_'+str(i)+'.csv')
    if i == 0:
        
        x = data['Points:0'].to_numpy()
        dom1 = np.where((x <= 2e-4) & (x >= 1e-4))[:][0]
        x1 = x[dom1]

        ars = np.zeros(shape=(nt, len(x1)))
    #else:
    #    break
    ars[i,:] = np.exp(data['mean_en'][dom1])*6.022e23
    time[i] = data['Time'][0]

xx, yy = np.meshgrid(x1, time)

plt.contourf(xx, yy, ars, 100)
plt.colorbar()
plt.show()
