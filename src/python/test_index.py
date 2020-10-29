import ctypes
import numpy as np

lib = ctypes.cdll.LoadLibrary('../cpp/libindex.so')
lib.c_get_mean_rank.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.c_get_mean_rank.restype = ctypes.POINTER(ctypes.c_float)
lib.c_free_mem.argtypes = [ctypes.POINTER(ctypes.c_float)]
lib.c_make_represent_matrix.argtypes = [ctypes.c_char_p,
                                        ctypes.c_char_p,
                                        ctypes.c_int,
                                        ctypes.POINTER(ctypes.c_int), 
                                        ctypes.c_int]


lib.c_convert_to_gene_ontology.argtypes = [ctypes.c_char_p,
                                          ctypes.c_char_p]

path = '/home/tung/RepresentData/GSE150430/matrix.hdf5'
path_out = '/home/tung/RepresentData/Transform/GSE150430'
lib.c_convert_to_gene_ontology(path.encode(), path_out.encode())

# string = b"//home//tung//.BioTBDataDev//Data//SingleCell//Study//TabulaMuris_droplet//main//matrix.hdf5"
# string_newhdf5 = b"/home/tung/new.hdf5"

# arr_cell = np.zeros(55657, dtype=np.int32)
# arr_cell[2000:] = 1


# pointer_arr_cell = arr_cell.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
# x = lib.c_make_represent_matrix(string,
#                                 string_newhdf5,
#                                 ctypes.c_int(55656),
#                                 pointer_arr_cell,
#                                 ctypes.c_int(2))


