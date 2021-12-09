from colvars_grid import *

histo = colvars_grid('DBC.dat')
import numpy as np
import matplotlib.pyplot as plt
plt.plot(np.transpose(*histo.meshgrid()),histo.data[0], label='DBC')
plt.legend()
plt.show()
