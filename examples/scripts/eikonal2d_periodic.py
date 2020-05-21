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
parser.add_argument('--target', type=float, nargs=2, default=[0.0, 0.0],
                    help='Coordinates of the target point.')
parser.add_argument('--no_plot', action='store_true',
                    help='Don\'t plot the result.')

args = parser.parse_args()

filename = '../data/eikonal2d_periodic.h5'

if os.path.isfile(filename):
    os.remove(filename)

f = h5.File(filename, "w")
data = np.ones((args.npts, args.npts))

i = int(np.floor((args.npts - 1) * (1 + args.target[0]) / 2.0))
j = int(np.floor((args.npts - 1) * (1 + args.target[1]) / 2.0))

data[i, j] = -1

f.create_dataset('cost_function', shape=data.shape, data=data)
f.close()

subprocess.run(["../../build/examples/eikonal2d_periodic"])

f = h5.File(filename, "r")
data = f['value_function'][()]
f.close()

t = np.linspace(-1, 1, data.shape[0])
x, y = np.meshgrid(t, t, indexing='ij')

if not args.no_plot:
    plt.pcolormesh(x, y, data, cmap='jet')
    plt.colorbar()
    plt.contour(x, y, data, levels=10, linestyles='dashed', colors='k')
    plt.axis('equal')
    plt.show()
