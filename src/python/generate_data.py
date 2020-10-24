from CS_data_processing import *
from CS_utils import *
import numpy as np
import traceback
import ctypes
import json
import os

lib = ctypes.cdll.LoadLibrary('../cpp/libindex.so')
lib.c_get_mean_rank.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
lib.c_get_mean_rank.restype = ctypes.POINTER(ctypes.c_float)
lib.c_free_mem.argtypes = [ctypes.POINTER(ctypes.c_float)]
lib.c_make_represent_matrix.argtypes = [ctypes.c_char_p,
                                        ctypes.c_char_p,
                                        ctypes.c_int,
                                        ctypes.POINTER(ctypes.c_int), 
                                        ctypes.c_int,
                                        ctypes.c_char_p]

def make_represent_matrix(path_hdf5,
                          path_hdf5_out,
                          number_cells,
                          cluster_id,
                          number_clusters,
                          study_name):
    number_cells = int(number_cells)
    number_clusters = int(number_clusters)
    cluster_id_int = np.array(cluster_id, dtype=np.int32)
    pointer_arr_cell = cluster_id_int.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    lib.c_make_represent_matrix(path_hdf5.encode(),
                                path_hdf5_out.encode(),
                                ctypes.c_int(number_cells),
                                pointer_arr_cell,
                                ctypes.c_int(number_clusters),
                                study_name.encode())

def index_study(study_dir, out_path_dir, study_id):
    species = get_study_species(study_dir)
    if species != "human":
        return False

    annotation = load_annotations(study_dir)

    if is_study_annotated(annotation):
        cell_based_cluster_label, cell_based_clusters = load_study_cell_type_cluster(study_dir)
    else:
        return False

    if len(cell_based_cluster_label) <= 4:
        return False

    create_path(out_path_dir)

    file = open(os.path.join(out_path_dir, "metadata.json"), "w")
    json.dump(cell_based_cluster_label, file)
    file.close()
    file = open(os.path.join(out_path_dir, "gene_list.json"), "w")
    json.dump(load_genes_list(study_dir), file)
    file.close()

    try:
        make_represent_matrix(os.path.join(study_dir, "main", "matrix.hdf5"),
                              out_path_dir,
                              len(cell_based_clusters),
                              cell_based_clusters,
                              len(cell_based_cluster_label),
                              study_id)
    except Exception as e:
        print(e)
        os.rmdir(out_path_dir)
        return False

    return True

def make_represent_data(path, path_out):
    create_path(path_out)
    studies_path = [name for name in os.listdir(path)]
    for study in studies_path:
        print(study)
        if study in ["GSE136831",
                    "MCA", "SCP498",
                    "GSE121611_cryobiopsy",
                    "ventotormo2018",
                    "MOCA_2M", "zhong2018"]:
            continue
        model_path = os.path.join(path_out, study)
        if os.path.exists(model_path):
            continue
        try:
            index_study(os.path.join(path, study), model_path, study)
        except Exception as e:
            print(e)
            traceback.print_exc()

make_represent_data('/home/tung/.BioTBDataDev/Data/SingleCell/Study', '/home/tung/RepresentData')


    

        






    
    