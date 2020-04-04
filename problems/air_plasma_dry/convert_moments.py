import numpy as np


diff_file = np.loadtxt('electron_diffusivity.txt')
mu_file = np.loadtxt('electron_mobility.txt')
N = 2.445692e19

f = 'electron_moments.txt'
with open(f,'w') as write_file:
    for i in range(0,len(diff_file)):
        write_file.write('{0:.6e} {1:.6e} {2:.6e} \n'.format(diff_file[i, 0], mu_file[i, 1]/N, diff_file[i, 1]/N))
    write_file.closed
