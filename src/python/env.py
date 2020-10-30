PATH_GO2SIZE = '../../data/go2size.json' 
PATH_GO2SIZE = '../../data/go2size.json'
PATH_GO_GENE = '../../data/go_gene.json'
PATH_GO2INDEX = '../../data/go2index.json'
PATH_ENS_HUMAN = '../../data/ens_human.json'
#PATH_GO_TREE = '../../data/goslim_generic.json'
PATH_GSE2INDEXGO = '../../data/gse2indexgo.json'
PATH_HUMAN_ENRICH = '../../data/human_enrich.json'
PATH_INDEXGO2SIZE = '../../data/indexgo2size.json'
PATH_GO_TREE = '../../data/gene_onto.txt'
PATH_MART = '../../data/mart_export.txt'

import json


def get_go2index():
    with open(PATH_GO2INDEX, 'r') as fi:
        return json.load(fi)

def get_index2go():
    ans = {}
    go2index = get_go2index()
    for go in go2index:
        ans[go2index[go]] = go
    return ans

def get_mart():
    ans = {}
    with open(PATH_MART, 'r') as fi:
        for line in fi:
            row = line.strip()
            data = row.split('\t')
            if len(data) != 2:
                continue
            if data[1] not in ans:
                ans[data[1]] = {'genes': []}
            ans[data[1]]['genes'].append(data[0])
    return ans




