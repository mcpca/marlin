#!/usr/bin/env python3

import numpy as np
import h5py as h5
import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--npts', type=int, default=101,
                    help='Number of gridpoints per dimension.')
parser.add_argument('--target', type=float, nargs=3, default=[0.0, 0.0, 0.0],
                    help='Coordinates of the target point.')

args = parser.parse_args()

filename = '../data/eikonal3d.h5'

if os.path.isfile(filename):
    os.remove(filename)

f = h5.File(filename, "w")
data = np.ones((args.npts, args.npts, args.npts))

i = int(np.floor((args.npts - 1) * (1 + args.target[0]) / 2.0))
j = int(np.floor((args.npts - 1) * (1 + args.target[1]) / 2.0))
k = int(np.floor((args.npts - 1) * (1 + args.target[2]) / 2.0))

data[i, j, k] = -1

f.create_dataset('cost_function', shape=data.shape, data=data)
f.close()

subprocess.run(["../../build/examples/eikonal3d"])

f = h5.File(filename, "r")

data = f['value_function'][()]
f.close()

t = np.linspace(-1, 1, data.shape[0])
x, y, z = np.meshgrid(t, t, t, indexing='ij')
w = np.sqrt((x - args.target[0])**2 + (y - args.target[1])**2 +
            (z - args.target[2])**2)

print('Maximum error:', np.amax(np.abs(w - data)))
