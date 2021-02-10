import numpy as np
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import json


with open('out.json', 'r') as fp:
    error = json.load(fp)

radius = [ 2, 3, 4, 5, 6 ]
xticks = radius + [ 7, 8 ]
markers = [ 'o', '*', '^', 's' ]

for key, marker in zip(error, markers):
    plt.plot(radius, error[key], marker=marker, label=key)

plt.ylim([0, 30])
plt.xticks(xticks)
plt.xlabel('Convolution kernel radius (pixels)', fontsize='xx-large')
plt.ylabel('Error (arbitrary units)', fontsize='xx-large')
plt.legend(title = 'Oversampling', fontsize='large', title_fontsize='x-large')
plt.savefig('error_v_kernelradius.png', dpi=300)

