import numpy as np
import h5py
import json
import os
import tensorflow as tf
from sklearn import preprocessing
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier, plot_tree
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from collections import Counter
import graphviz 
from load_generate_data import *



#profile = load_study(study_path)

def get_highly_genes_in_study(study_path):
    profile = load_study(study_path)
    standardized_profile = profile
    cluster_profile = standardized_profile
    cluster_centroid = np.array(np.mean(cluster_profile, axis = 0).tolist()[0], dtype=np.float32)
    highly_genes_indices = np.transpose(np.argwhere(cluster_centroid > np.mean(cluster_centroid))).tolist()[0]
    map_id_index = load_genes_id2index_mapping(study_path)
    highly_in_study = [map_id_index['index'][index] for index in highly_genes_indices]
    return highly_in_study

def get_highly_genes_in_list_study(list_study_path):
    highly_genes = get_highly_genes_in_study(list_study_path[0])

    # Get intersection of highly genes express in all study
    for study_path in list_study_path[1:]:
        highly_in_study = get_highly_genes_in_study(study_path)
        highly_genes = list(set(highly_genes) & set(highly_in_study))
        
    return highly_genes

def map_highly_genes_id_to_study_col(highly_genes, study_path):
    map_id_index = load_genes_id2index_mapping(study_path)
    index = np.where([map_id_index['id'][id] for id in highly_genes] != None)
    return index


def get_rank_matrix(sparse_matrix):
    for cell_id in range(sparse_matrix.shape[0]):
        row_index = sparse_matrix[cell_id, :].indices
        data_row = sparse_matrix[cell_id, :].data
        sparse_matrix[cell_id, row_index] = np.argsort(data_row) + 1
    return sparse_matrix


def get_data_joined(path):

    list_study_path =   ['/home/tung/RepresentData/Transform/GSE110686',
                        '/home/tung/RepresentData/Transform/GSE124310',
                        '/home/tung/RepresentData/Transform/GSE150430']
                    # ['/home/tung/RepresentData/Transform/GSE123904_human',
                    # '/home/tung/RepresentData/Transform/GSE135922_all',
                    # '/home/tung/RepresentData/Transform/GSE140228_10X',
                    # '/home/tung/RepresentData/Transform/GSE150430']
    data = None
    have_data = False
    barcodes = []
    for study in list_study_path:
        mt, barcode_study = load_generate_data(study)
        print(mt.shape)
        if have_data == False:
            have_data = True
            data = mt
        else:
            data = np.concatenate((data, mt), axis=0)
        barcodes.extend(barcode_study)

        


    # with open(os.path.join(path, "data.json")) as fi:
    #     data = json.load(fi)
    # with open(os.path.join(path, "meta.json")) as fi:
    #     barcodes = json.load(fi)
    
    barcode_origin = np.array([('_').join(x.split('_')[1:]) for x in barcodes])
    barcodes = [x.split('_')[-1] for x in barcodes]
    le = preprocessing.LabelEncoder()                                                                                                                                         
    le.fit(barcodes)                                                                                                                                                          
    barcodes = le.transform(barcodes)
    print(barcodes)

    # kmeans = KMeans(n_clusters=50, random_state=0).fit(data)
    # group = np.array(kmeans.labels_)

    # for i in range(max(group) + 1):
    #     list_index = np.where(group == i)
    #     print(Counter(barcode_origin[list_index]))
    #     print("########################################")
    #     print("########################################")
    #     print("########################################")
    #     print("########################################")

    # print()
    
    data = np.array(data)

    value = sum(data > 0.5, 1)
    print(len(value))
    number_cells = data.shape[0]
    print(number_cells / 2)
    list_good = np.where(value < (number_cells / 2))[0]
    index_good = np.array(range(data.shape[1]))[list_good]


    clf = tree.DecisionTreeClassifier(max_depth=10, max_leaf_nodes=50)
    clf = clf.fit(data, barcodes)
    print(clf)
    print(clf.score(data, barcodes))
    plot_tree(clf, filled=True)
    plt.show()

    dot_data = tree.export_graphviz(clf, out_file=None) 
    graph = graphviz.Source(dot_data) 
    graph.render("iris") 




    # data = np.expand_dims(data, -1)

    # model = tf.keras.Sequential([
    #     tf.keras.layers.Flatten(input_shape=(20,1)),
    #     tf.keras.layers.Dense(64, activation='elu'),
    #     tf.keras.layers.Dense(128, activation='elu'),
    #     tf.keras.layers.Dense(64, activation='elu'),
    #     tf.keras.layers.Dense(50),
    #     tf.keras.layers.Softmax()
    # ])
    # model.summary()
    # model.compile(optimizer='adam',``
    #             loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
    #             metrics=['accuracy'])
    # model.fit(data, barcodes, epochs=100)




#get_data_joined('/home/tung/RepresentData/Transform')


    

    



