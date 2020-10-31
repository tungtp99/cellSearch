import json
import matplotlib.pyplot as plt
from sklearn import preprocessing
import os
import numpy as np


def get_name(study_name, lable):
    if lable.find("BTC000220") != -1:
        return(study_name + lable)
    return ''

def plot_tsne_liger(path):
    map_id_celltype_name = {'U': 'Unassigned'}
    def get_color(name):
        map_id_celltype_name[name.split('_')[1]] = name.split('_')[1]
        return name.split('_')[1]
        if name.find('BTC') != -1:
            id = name.split('BTC')[1].split(']')[0]
            map_id_celltype_name[id] = name.split('_')[-1]
            return id
        return 'U'

    with open(os.path.join(path,'meta.json'), 'r') as fi:
        meta = json.load(fi)

    with open(os.path.join(path,'tsne.json'), 'r') as fi:
        tsne = json.load(fi)

    x = [x[0] for x in tsne]
    y = [x[1] for x in tsne]

    barcodes = [get_color(x) for x in meta]
    origin_barcodes = np.array(barcodes)
    # print(barcodes)
    # print(set(barcodes))


    def onpick3(event):
        ind = event.ind
        print('onpick3 scatter:', [map_id_celltype_name[i] for i in list(set(origin_barcodes[ind]))])


    le = preprocessing.LabelEncoder()
    le.fit(barcodes)
    barcodes = le.transform(barcodes)
    print(barcodes)


    fig, ax = plt.subplots()
    ax.scatter(x, y, c=barcodes, cmap=plt.get_cmap('tab20c'), picker=True)
    fig.canvas.mpl_connect('pick_event', onpick3)
    plt.show()


#plot_tsne_liger('/home/tung/RepresentData/Transform')