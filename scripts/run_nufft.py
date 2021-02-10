import numpy as np
from skimage.data import shepp_logan_phantom
import subprocess
import json


data = shepp_logan_phantom()
norm = np.sum(data * data)
ncol = 400
nrow = 400

def calc_error(output):
    scale = np.sum(output * data) / norm
    return np.linalg.norm(output - data * scale)

if __name__ == '__main__':

    over_samples = [ '1.5', '2', '3', '4' ]
    kern_radius = [ '2', '3', '4', '5', '6' ]

    prog = ['./nufft' ]
    error = { key : [] for key in over_samples }

    for o in over_samples:
        error[o] = []
        for r in kern_radius:
            prog_args = prog + [o, r]
            print('running ...', str(prog_args))
            p = subprocess.Popen(prog_args)
            rc = p.wait()
            if not rc == 0:
                print('something went wrong', rc);
                exit(1)
            output = np.fromfile('output.bin', dtype=np.float32)
            output = output.reshape(nrow, ncol)
            e = calc_error(output);
            print('error ', e)
            error[o].append(e)
            fname = 'shepp_%s_%s.npy' % (o, r)
            np.save(fname, output)

    with open('error.json', 'w') as fp:
        json.dump(error, fp)
