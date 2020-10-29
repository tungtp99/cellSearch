from env import *
import json


def parese_go_tree():
    with open(PATH_GO_TREE, 'r') as fi:
        tree_structer = json.load(fi)

    print(tree_structer['graphs'])

parese_go_tree()
    


