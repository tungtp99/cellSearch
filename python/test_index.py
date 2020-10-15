import ctypes

lib = ctypes.cdll.LoadLibrary('../libindex.so')

lib.c_get_mean_rank.argtypes = [ctypes.c_char_p]
lib.c_get_mean_rank.restype = ctypes.POINTER(ctypes.c_float)

lib.c_free_mem.argtypes = [ctypes.POINTER(ctypes.c_float)]


string = b"//home//tung//.BioTBDataDev//Data//SingleCell//Study//TabulaMuris_droplet//main//matrix.hdf5"

x = lib.c_get_mean_rank(string)
lib.c_free_mem(x)

