import json
import os
import pathlib

import h5py
import numpy as np

def load_sparse_data(study_path):
    # Read study metadata
    f = h5py.File(os.path.join(study_path, 'main', 'matrix.hdf5'), 'r')
    
    # The naming of each gene and each cell
    cells_code = f['bioturing']['barcodes']
    genes_code = f['bioturing']['features']

    # Total genes and cells in the study
    genes_count = len(genes_code)
    cells_count = len(cells_code)
    
    # Data matrix
    data = f['bioturing']['data']
    data = data[0:len(data)]

    indices = f['bioturing']['indices'] # Gene index - row
    indices = indices[0:len(indices)]
    
    indptr = f['bioturing']['indptr']
    
    return data, indices, indptr, cells_count, genes_count

def load_sparse_data_row_col(study_path):
    f = h5py.File(os.path.join(study_path, 'main', 'matrix.hdf5'), 'r')

    # Data matrix
    data = f['bioturing']['data']
    data = data[0:len(data)]

    indices = f['bioturing']['indices'] # Gene index - row
    indices = indices[0:len(indices)]
    
    indptr = f['bioturing']['indptr']

    ans = []
    pos = 0
    for cell_index, _ in enumerate(indptr):
        for i, row_index in indices[indptr[cell_index]:indptr[cell_index + 1]]:
            ans.append([cell_index, row_index, data[pos + i]])
        pos += indptr[cell_index + 1] - indptr[cell_index]
    return ans

def load_gene_sparse_data(study_path):
    # Read study metadata
    f = h5py.File(os.path.join(study_path, 'main', 'matrix.hdf5'), 'r')
    
    # The naming of each gene and each cell
    cells_code = f['bioturing']['barcodes']
    genes_code = f['bioturing']['features']

    # Total genes and cells in the study
    genes_count = len(genes_code)
    cells_count = len(cells_code)
    
    # Data matrix
    data = f['countsT']['data']
    data = data[0:len(data)]

    indices = f['countsT']['indices'] # Gene index - row
    indices = indices[0:len(indices)]
    
    indptr = f['countsT']['indptr']
    
    return data, indices, indptr, cells_count, genes_count
    

def load_study(study_path):
    all_cells = []
    
    # Read study metadata
    f = h5py.File(os.path.join(study_path, 'main', 'matrix.hdf5'), 'r')
    
    # The naming of each gene and each cell
    cells_code = f['bioturing']['barcodes']
    genes_code = f['bioturing']['features']

    # Total genes and cells in the study
    genes_count = len(genes_code)
    print(genes_count)
    cells_count = len(cells_code)
    
    # Data matrix
    data = f['bioturing']['data']
    data = data[0:len(data)]

    indices = f['bioturing']['indices'] # Gene index - row
    indices = indices[0:len(indices)]
    
    indptr = f['bioturing']['indptr']

    from scipy.sparse import csc_matrix
    return csc_matrix((data, indices, indptr), shape=(genes_count, cells_count)).transpose()

def load_cells_at_index(study_path, cell_indices):
    all_cells = []
    
    # Read study metadata
    f = h5py.File(os.path.join(study_path, 'main', 'matrix.hdf5'), 'r')
    
    # The naming of each gene and each cell
    cells_code = f['bioturing']['barcodes']
    genes_code = f['bioturing']['features']

    # Total genes and cells in the study
    genes_count = len(genes_code)
    cells_count = len(cells_code)
    
    # Data matrix
    data = f['bioturing']['data']
    data = data[0:len(data)]

    indices = f['bioturing']['indices'] # Gene index - row
    indices = indices[0:len(indices)]
    
    indptr = f['bioturing']['indptr']
    
    cells = []
    ranges = list(map(lambda x: (x, x + 1), cell_indices))
    ranges = list(map(lambda x: (indptr[x[0]], indptr[x[1]]), ranges))
    for _range in ranges:
        cell = [0] * genes_count
        vals = data[_range[0] : _range[1]]
        idx = indices[_range[0] : _range[1]]
        for i in range(0, len(idx)):
            cell[idx[i]] = vals[i]

        cells.append(cell)

    return np.array(cells)

def load_genes_at_index(gene_indices, data, indices, indptr, cells_count, genes_count):
    genes = []
    ranges = list(map(lambda x: (x, x + 1), gene_indices))
    ranges = list(map(lambda x: (indptr[x[0]], indptr[x[1]]), ranges))
    for _range in ranges:
        gene = [0] * cells_count
        vals = data[_range[0] : _range[1]]
        idx = indices[_range[0] : _range[1]]
        for i in range(0, len(idx)):
            gene[idx[i]] = vals[i]

        genes.append(gene)

    return np.array(genes)

def load_cells_at_index(cell_indices, data, indices, indptr, cells_count, genes_count):
    cells = []
    ranges = list(map(lambda x: (x, x + 1), cell_indices))
    ranges = list(map(lambda x: (indptr[x[0]], indptr[x[1]]), ranges))
    for _range in ranges:
        cell = [0] * genes_count
        vals = data[_range[0] : _range[1]]
        idx = indices[_range[0] : _range[1]]
        for i in range(0, len(idx)):
            cell[idx[i]] = vals[i]

        cells.append(cell)

    return np.array(cells)

def load_study_graph_based_cluster(study_path):
    data = load_annotations(study_path)
    result = {}
    return data["graph_based"]["clusters"]

def load_study_cell_type_cluster(study_path):
    annotation = load_annotations(study_path)
    for key in annotation.keys():
        if "bioturingcelltype" in (annotation[key]["name"]).lower().replace('-', '').replace(' ', ''):
            clusterName = annotation[key]["clusterName"]
            cluster = annotation[key]["clusters"]
            return clusterName, cluster
    return None, None

def is_study_annotated(annotation):
    for key in annotation.keys():
        if "bioturing" in (annotation[key]["name"]).lower():
            return True
    return False

def get_study_cells_in_types(annotation, cells_indices = None):
    for key in annotation.keys():
        if "bioturingcelltype" in (annotation[key]["name"]).lower().replace('-', '').replace(' ', ''):
            
            cell_types = list(map(lambda x: annotation[key]["clusterName"][x].strip(), annotation[key]["clusters"]))
            if cells_indices is not None:
                cell_types = list(map(lambda x: x.strip(), np.array(cell_types)[cells_indices]))
            return cell_types
    return None

def get_study_info(study_path):
    study_info = open(os.path.join(study_path, "run_info.json")).read()
    json_data = json.loads(study_info)
    
    return json_data

def get_study_species(study_path):
    study_info = get_study_info(study_path)
    species = study_info.get("index_type")
    return species.lower()

def get_study_tags(study_path):
    study_info = get_study_info(study_path)
    tags = study_info.get("tag")
    if tags is None:
        return None
    return list(map(lambda x: str(x).lower(), tags))

def load_annotations(study_path):
    f = open(os.path.join(study_path, "main", "metadata", "metalist.json"))
    metalist_data = json.loads(f.read())
    
    for cluster_id in metalist_data.keys():
        f = open(os.path.join(study_path, "main", "metadata", cluster_id + ".json"))
        metalist_data[cluster_id]["clusters"] = json.loads(f.read())
    
    return metalist_data

def load_cluster_result(study_path):
    data = load_annotations(study_path)
    result = {"graph": {}}
    result["graph"]["clusters"] = data["graph_based"]["clusters"]

    return result

def load_genes_id2index_mapping(study_path):
    # Perform the mapping
    genes_list = load_genes_list(study_path)
    mapping = {}
    mapping["id"] = {gene_id : index for index, gene_id in enumerate(genes_list)}
    mapping["index"] = {index : gene_id for index, gene_id in enumerate(genes_list)}
    
    return mapping

def load_genes_list(study_path):
    f = h5py.File(os.path.join(study_path, 'main', 'matrix.hdf5'), 'r')
    genes_list = []
    try:
        genes_list = list(map(lambda x: x.decode("utf-8"), list(f['bioturing']['features'])))
    except Exception as e:
        genes_list = list(f['bioturing']['features'])
    
    return genes_list
    
def create_path(path):
    # Create paths leading to the file
    os.makedirs(path, exist_ok=True)

