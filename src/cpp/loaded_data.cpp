#include "loaded_data.h"
#include "hdf5io.h"


string pathAppend(const string& p1, const string& p2) {
  char sep = '/';
  string tmp = p1;

#ifdef _WIN32
  sep = '\\';
#endif

  if (p1[p1.length()] != sep) { // Need to add a
    tmp += sep;                 // path separator
    return(tmp + p2);
  }
  return(p1 + p2);
}

struct loaded_data_t *load_data(const string *path,
                                int use_reduce)
{
	auto res = (struct loaded_data_t *) malloc(sizeof(struct loaded_data_t));
	File file(*path);
	Group group;
	group = file.getGroup("bioturing");
  fprintf(stdout, "Loading data\n");
  res->shape = extract1DArray<int>(group, "shape");
  res->data = extract1DArray<float>(group, "data");
  res->barcodes = extract1DArray(group, "features");
  res->gene_names = extract1DArray(group, "barcodes");
  res->indices = extract1DArray<int>(group, "indices");
	res->indptr = extract1DArray<long unsigned>(group, "indptr");
  return res;
}

struct loaded_data_t *load_data(const string *path,
                                const string *group_hdf5,
                                int use_reduce)
{
	auto res = (struct loaded_data_t *) malloc(sizeof(struct loaded_data_t));
	File file(*path);
	Group group;
	group = file.getGroup(*group_hdf5);
  fprintf(stdout, "Loading data\n");
  res->shape = extract1DArray<int>(group, "shape");
  res->data = extract1DArray<float>(group, "data");
  res->barcodes = extract1DArray(group, "features");
  res->gene_names = extract1DArray(group, "barcodes");
  res->indices = extract1DArray<int>(group, "indices");
	res->indptr = extract1DArray<long unsigned>(group, "indptr");
  return res;
}

void write_hdf5(const std::string* path,
                vvvp matrix,
                const std::string* study_name)
{
  std::vector<float> data;
  std::vector<int> indices;
  std::vector<long unsigned> indptr = {0};
  std::vector<std::string> barcodes;
  
  std::ifstream i;
  i.open(pathAppend(*path, "metadata.json"));
  json cell_type_json;
  i >> cell_type_json;
  i.close();
  i.open(pathAppend(*path, "gene_list.json"));
  json list_gene_json;
  i >> list_gene_json;
  i.close();

  std::vector<std::string> cell_type = cell_type_json.get<std::vector<std::string>>();
  std::vector<std::string> list_gene = list_gene_json.get<std::vector<std::string>>();
  std::cout << "Number of cluster is " << matrix.size() << "\n";
  
  int pos = 0;
  int number_cells = 0;
  for (int cluster; cluster < matrix.size(); ++cluster) {
    for (int c = 0; c < matrix[cluster].size(); ++c) {
      std::sort(matrix[cluster][c].begin(),
                matrix[cluster][c].end(),
                [](std::pair<float, int> a, std::pair<float, int> b) {
                  return a.second < b.second; 
                });
      for (int g = 0; g < matrix[cluster][c].size(); ++g) {
        indices.push_back(matrix[cluster][c][g].second);
        data.push_back(matrix[cluster][c][g].second);
        pos++;
      }
      indptr.push_back(pos);
      barcodes.push_back(std::to_string(c) + "_"
                        + *study_name + "_"
                        + cell_type[cluster]);
    }
    number_cells += matrix[cluster].size();
  }

  std::vector<int> shape = {int(list_gene.size()), number_cells};

  write_matrix_hdf5("/", *path, &data, &shape,
                    &indptr, &indices, &barcodes,
                    &list_gene);
  
}


void write_matrix_hdf5(std::string group_hdf5,
                      std::string path,
                      std::vector<float> *data,
                      std::vector<int> *shape,
                      std::vector<long unsigned> *indptr,
                      std::vector<int> *indices,
                      std::vector<std::string> *barcodes,
                      std::vector<std::string> *gene_names) 
{
  File file(pathAppend(path, "matrix.hdf5"), File::ReadWrite | File::Create | File::Truncate);
  Group group = file.getGroup(group_hdf5);
  write1DArray(group, "shape", *shape);
  write1DArray(group, "indptr", *indptr);
  write1DArray(group, "indices", *indices);
  write1DArray<float>(group, "data", *data);
  write1DArray(group, "barcodes", *barcodes);
  write1DArray(group, "features", *gene_names);
}

std::unordered_map<int, int> transform_features_to_go(std::vector<std::string> features,
                              json &gse2indexgo,
                              std::unordered_map<std::string, std::vector<int>> &gse2indexgo_stl,
                              std::unordered_map<std::string, bool> &bad_features)
{
    std::unordered_map<int, int> go_value;
    for (int i = 0; i < features.size(); ++i) {
        if (gse2indexgo_stl.find(features[i]) == gse2indexgo_stl.end()) {
          if (bad_features.find(features[i]) != bad_features.end()) {
            continue;
          }
          if (gse2indexgo.find(features[i]) == gse2indexgo.end()) {
            bad_features[features[i]] = true;
            continue;
          }
          std::vector<int> go = gse2indexgo[features[i]].get<std::vector<int>>();
          gse2indexgo_stl[features[i]] = go;
        }
        
        if (gse2indexgo_stl.find(features[i]) == gse2indexgo_stl.end())
          continue;
        
        std::vector<int> go = gse2indexgo_stl[features[i]];
        for (int j = 0; j < go.size(); ++ j) {
          go_value[go[j]]++;
        }
    }
    return go_value;
}

void destroy_loaded_data(struct loaded_data_t *loaded_data)
{
  delete loaded_data->data;
  delete loaded_data->shape;
  delete loaded_data->indptr;
  delete loaded_data->indices;
  delete loaded_data->barcodes;
  delete loaded_data->gene_names;
}


