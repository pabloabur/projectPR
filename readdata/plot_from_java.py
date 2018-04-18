import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt('force_RC.txt', unpack=True)
plt.plot(x, y, '-')
plt.ylabel('Forca')
plt.ylabel('Tempo (ms)')
plt.show()
