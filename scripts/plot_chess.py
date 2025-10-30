import sys
import numpy as np
import matplotlib.pyplot as plt
import tifffile 

if __name__ == '__main__':
    if len(sys.argv) < 2:
        exit(1)

    data = tifffile.imread(sys.argv[1]);
    nslcs = data.shape[0]
    idxs = [i for i in range(nslcs) if data[i].max() > 0]
    idxs = [idxs[0], nslcs//2, idxs[-1]]
    fig, axs = plt.subplots(ncols = len(idxs), figsize = [10, 5])
    for i, j in enumerate(idxs):
        axs[i].imshow(data[j], cmap='gray')
    plt.show()
        
