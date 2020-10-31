from load_generate_data import load_generate_data
import matplotlib.pyplot as plt
import networkx as nx
from env import *
import json
import csv


def get_id(link):
    return link.split('/')[-1].replace('_', ':')

def parese_go_tree():
    with open(PATH_GO_TREE, 'r') as fi:
        tree_structer = json.load(fi)
    name_go = {}
    link = {}
    
    for go in tree_structer:
        u = go['id']
        if 'name' in go:
            name_go[u] = go['name']
        else:
            name_go[u] = u
        link[u] = {}
        if 'child' not in go:
            continue
        if u not in link:
            link[u] = {}
        for v in go['child']:
            link[u][v] = 'is_a'
        
    return link, name_go

def draw_graph(path_study, path_graph, cell_index):
    def get_name(id, name_go):
        if id not in name_go:
            return id
        return name_go[id].split('/')[-1] + id

    mt, barcodes =  load_generate_data(path_study, )
    print("get data")
    print(barcodes[cell_index])
    link, name_go = parese_go_tree()
    print("get gene ontology structer")

    cell_data = mt[cell_index,:]
    print(mt.shape)
    index2go = get_index2go()
    list_go = set()
    for index, value in enumerate(cell_data):
        if value >= 0.5:
            go = index2go[index]
            list_go.add(go)
    
    graph = [['id_from', 'id_to', 'data']]
    for go in list_go:
        for v in list_go:
            if go in link and v in link[go]:
                graph.append([get_name(go, name_go) , get_name(v, name_go), link[go][v]])

    with open(path_graph, 'w') as fi:
        writer = csv.writer(fi)
        writer.writerows(graph)            