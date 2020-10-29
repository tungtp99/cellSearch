import json
import os

def get_geneId(path):
    with open(os.path.join(path, "gene_list.json")) as fi:
        data = json.load(fi)
    return data

def intersection(list1, list2):
    return list(set(list1) & set(list2))
    
def get_list_geneId(path_studies, path_study_names):
    list_studies = [i.split('.txt')[0] for i in os.listdir(path_study_names)]

    map_list_gene = {}
    for studyId in list_studies:
        path_to_study = os.path.join(path_studies, studyId)
        map_list_gene[studyId] = get_geneId(path_to_study)
    for x in map_list_gene:
        for y in map_list_gene:
            print(len(intersection(map_list_gene[x], map_list_gene[y])), x, y)

get_list_geneId('/home/tung/RepresentData', '/home/tung/220')



