import ctypes
import numpy as np

lib = ctypes.cdll.LoadLibrary('../cpp/libindex.so')
lib.c_get_mean_rank.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.c_get_mean_rank.restype = ctypes.POINTER(ctypes.c_float)
lib.c_free_mem.argtypes = [ctypes.POINTER(ctypes.c_float)]
string = b"//home//tung//.BioTBDataDev//Data//SingleCell//Study//TabulaMuris_droplet//main//matrix.hdf5"

arr_cell = np.array([i for i in range(50000)], dtype = np.int32)
pointer_arr_cell = arr_cell.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

x = lib.c_get_mean_rank(string, pointer_arr_cell, ctypes.c_int(len(arr_cell)))
tmp = np.array(x[0:23433])
lib.c_free_mem(x)

import matplotlib.pyplot as plt 

fig, ax = plt.subplots()
plt.ylim(0, 23433)
ax.plot(np.sort(tmp))
plt.show()

