# Simple script to convert *.h5 files output by Mathematica to *.npy files, which are
# faster and easier to lead by Python/Numpy/Cython.

import h5py
import numpy
import os

all_file_names = os.listdir('.')

for file_name in all_file_names:
    if os.path.splitext(file_name)[1] == '.h5':
        file = h5py.File(file_name, 'r')
        data = file.get(file.keys()[0]).value
        numpy.save(os.path.splitext(file_name)[0] + '.npy', data)
