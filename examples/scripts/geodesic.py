#!/usr/bin/env python3

import numpy as np
import h5py as h5
import argparse
import os
import subprocess
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--npts', type=int, default=201,
                    help='Number of gridpoints per dimension.')

args = parser.parse_args()

filename = '../data/geodesic.h5'

if os.path.isfile(filename):
    os.remove(filename)

f = h5.File(filename, "w")
data = np.ones((args.npts, args.npts))

data[args.npts//2, args.npts//2] = -1

f.create_dataset('cost_function', shape=data.shape, data=data)
f.close()

subprocess.run(["../../build/examples/geodesic"])

f = h5.File(filename, "r")

data = f['value_function'][()]
f.close()

t = np.linspace(-0.5, 0.5, data.shape[0])
x, y = np.meshgrid(t, t, indexing='ij')

plt.pcolormesh(x, y, data, cmap='jet')
plt.colorbar()
plt.contour(x, y, data, levels=25, linestyles='dashed', colors='k')
plt.axis('equal')
plt.show()
