import os
import joblib
import numpy as np
DEBUG = False


os.environ['API_DATA'] = "/home/tung/Thesis/Index"
os.environ['STUDIES_DATA'] = "/mnt/hdd2/api_server"


### ['clusters', 'clusters_meta', 'species', 'id2index_mapping' 'mean_rank']
### Structure of cluster ['centroid', 'genes_indices', 'cluster_size', 'cluster_name']
def load_all():
  combie_cell_type = {}
  __models_path = os.path.join(os.environ.get('API_DATA'),'cell_search_models')
  model_folders = [name for name in os.listdir(__models_path)]

  for i, model_id in enumerate(model_folders):
      model_path = os.path.join(__models_path, model_id, "model_{}.pkl".format(model_id))
      model = joblib.load(model_path)
      index2id = model['id2index_mapping']['index']
      index_mapping = np.zeros(len(index2id.items()), dtype=object)
      for i, id in index2id.items():
          index_mapping[i] = id
      if model['species'] == 'mouse':
          continue

      clusters_meta = model['clusters_meta']
      for cluster in clusters_meta:
          name = cluster['cluster_name']
          if name.find('[BTC') != -1:
              name = name.split('[')[1][:-1]
              combie_cell_type[name + '_'  + model_id] = {
                  "index2id": index2id,
                 "mean_rank": cluster['mean_rank'],
              }

  return combie_cell_type


'''
    combie_cell_type have structer of set 
    { 
        id = BTC..._<study name> : {
            index2id
            mean_rank
        }
    }

    Function get all id of index2id and return a new id2index of all
    data
'''
def load_all_geneId(combie_cell_type):
    id2index = {}
    count = 0

    for name in combie_cell_type:
        data = combie_cell_type[name]
        for index in data['index2id']:
            if data['index2id'][index] not in id2index:
                id2index[data['index2id'][index]] = count
                count += 1
    return id2index

''' 
    Function get a list of value with of geneId 
    and map to a list if index in all geneId
'''
def map_to_all_geneId(data, geneIds, all_geneId_map):
    ans = np.zeros(len(all_geneId_map))
    ans[np.array([all_geneId_map[i] for i in geneIds])] = data
    return ans






          