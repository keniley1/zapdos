import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def func(x, a, b, c):
    return(a * np.exp(-b*x) + c)

#data = pd.read_csv('temperature_along_centerline.csv')
#data = pd.read_csv('temperature_along_centerline_finer.csv')
data = pd.read_csv('centerline_data_1100.csv')

x = data['arc_length']
t = data['T']

# Save data to file
#data[['arc_length', 'T']].to_csv('test.csv')
header = ['arc_length', 'T']
data.to_csv("test_1100.csv", columns = header, header=False, index=False)
exit()

#popt, pcov = curve_fit(func, x, t)

plt.semilogy(x, t, 'k*')
#plt.plot(x, func(x, *popt), 'r')
plt.show()


