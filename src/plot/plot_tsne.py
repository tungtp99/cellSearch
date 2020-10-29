name_SDY = ["Unassigned","effector memory CD4 T cell [BTC000232]","central memory CD4 T cell [BTC000245]","memory CD8 T cell [BTC000150]","CD8 T cell [BTC000221]","NK cell [BTC000015]","Treg cell [BTC000222]","B cell [BTC000020]","CD4 T cell [BTC000220]","naive B cell [BTC000209]","plasmacytoid dendritic cell [BTC000025]","epithelial cell [BTC000114]","macrophage [BTC000011]","conventional dendritic cell [BTC000026]"]
name_ST56 = ["Unassigned","pericyte [BTC000038]","melanocyte cell [BTC000152]","fibroblast [BTC000049]","CD8 T cell [BTC000221]","endothelial cell [BTC000036]","retinal pigment epithelial cell [BTC000146]","macrophage [BTC000011]","Schwann cell [BTC000057]","NK cell [BTC000015]","mast cell [BTC000021]","CD4 T cell [BTC000220]","B cell [BTC000020]","natural killer T cell [BTC000065]"]


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
        #return name.split('_')[1]
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




        

plot_tsne_liger('/home/tung/RepresentData/Transform')



# silhouette

        


