from sklearn.decomposition import PCA
from bioinfokit.visuz import cluster
from scipy.sparse import csc_matrix
from sklearn.manifold import TSNE
from sklearn import preprocessing
import matplotlib.pyplot as plt
from collections import Counter
from scipy import sparse
import numpy as np
from env import *


import anndata
import json
import h5py
import os

def load_generate_data(study_dir):
    h5_data = h5py.File(os.path.join(study_dir, 'matrix.hdf5'), 'r')
    barcodes = h5_data['barcodes']
    indices = h5_data['indices']
    indptr = h5_data['indptr']
    data = h5_data['data']
    sparse_matrix = csc_matrix((data, indices, indptr))
    return sparse_matrix.transpose().toarray(), barcodes


#load_generate_data('/home/tung/RepresentData/GSE124310')

def get_list_gene_repersent(study_path, group):
    f = h5py.File(os.path.join(study_path, 'matrix.hdf5'), 'r')
    genes_list = []
    try:
        genes_list = list(map(lambda x: x.decode("utf-8"), list(f[group]['features'])))
    except Exception as e:
        genes_list = list(f[group]['features'])
    return genes_list

def share_genes(list_study):
    genes_total = get_list_gene_repersent(os.path.join('/home/tung/RepresentData', list_study[0]), '/')
    for study in list_study:
        genes = get_list_gene_repersent(os.path.join('/home/tung/RepresentData', study), '/')
        genes_total = list(set(genes).intersection(set(genes_total)))
    print(len(genes_total))
### 'GSE124310', 'GSE140228_smartseq2', 'GSE134809', 'GSE128169', 'GSE99254', 

#share_genes(['SDY997_landscape', 'GSE148837', '56535e196ac6-5c3e-80e7-f408e40e246a', 'GSE124310'])


def transform_features_to_go(list_features, gse2indexgo):
    go_value = {}
    for gene in list_features:
        if gene not in gse2indexgo:
            gse2indexgo[gene] = []
    a = np.array([go_index for gene in list_features for go_index in gse2indexgo[gene]])
    print(len(a))
    return Counter(a)
    for gene in list_features:
        #print(gene)
        if gene not in gse2indexgo:
            continue
        for go_index in gse2indexgo[gene]:
            if go_index not in go_value:
                go_value[go_index] = 0
            go_value[go_index] += 1
    return go_value

def load_generate_data_go(study_dir):
    study_dir = '/home/tung/RepresentData/AN1801'  
    h5_data = h5py.File(os.path.join(study_dir, 'matrix.hdf5'), 'r')
    barcodes = [x.split('_')[-1] for x in h5_data['barcodes']]
    features = np.array(h5_data['features'])
    indices = np.array(h5_data['indices'])
    indptr = h5_data['indptr']
    shape = h5_data['shape']
    data = h5_data['data']
    with open(PATH_GO2SIZE, 'r') as fi:
        go2size = json.load(fi)
    with open(PATH_INDEXGO2SIZE, 'r') as fi:
        indexgo2size = json.load(fi)

    with open(PATH_GSE2INDEXGO, 'r') as fi:
        gse2indexgo = json.load(fi)


    I = []
    J = []
    V = []
    num_go = 0
    print(indexgo2size)


    for i in range(shape[1]):
        print(i)
        
        list_features = features[indices[indptr[i]:indptr[i + 1]]]
        #for index in range(indptr[i], indptr[i + 1]):
        #    list_features[index - indptr[i]] = features[indices[index]]
        go_value = transform_features_to_go(list_features, gse2indexgo)
        for go in go_value:
            if go in go_value:
                I.append(i)
                J.append(go)
                V.append(go_value[go] / indexgo2size[str(go)])
                num_go = max(num_go, go)

    le = preprocessing.LabelEncoder() 
    le.fit(barcodes)
    barcodes = le.transform(barcodes)

    sparse_matrix = sparse.coo_matrix((V,(I,J)),shape=(shape[1],num_go + 1)).toarray()

    print(sparse_matrix)
    pca = PCA(n_components=20)
    principalComponents = pca.fit_transform(sparse_matrix)
    print(pca.explained_variance_ratio_)

    Y = TSNE().fit_transform(principalComponents)

    fig, ax = plt.subplots()
    ax.scatter(Y[:,0], Y[:,1], c=barcodes)
    plt.show()

    #load_generate_data_go('/home/tung/RepresentData/AN1801')

    b = csc_matrix((data, indices, indptr), shape=(shape[1], shape[0])).toarray()