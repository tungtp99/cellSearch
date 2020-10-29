#ifndef __LOADED_DATA__
#define __LOADED_DATA__
#include <vector>
#include <string>
#include <fstream>
#include "json.hpp"


#define PATH_HUMAN_ENRICH "../../data/human_enrich.json"
#define PATH_GSE2INDEXGO "../../data/gse2indexgo.json"
#define PATH_INDEXGO2SIZE "../../data/indexgo2size.json"
#define PATH_GO2SIZE "../../data/go2size.json"
#define PATH_GO_GENE "../../data/go_gene.json"
#define PATH_ENS_HUMAN "../../data/ens_human.json"
#define PATH_GO2INDEX "../../data/go2index.json"



typedef std::vector<std::pair<float, int>> vp;
typedef std::vector<std::vector<std::pair<float, int>>> vvp;
typedef std::vector<std::vector<std::vector<std::pair<float, int>>>> vvvp;

using namespace std;
using json = nlohmann::json;

struct loaded_data_t {
  json *annotation;
  vector<int> *shape;
  vector<float> *data;
  vector<int> *indices;
  vector<string> *barcodes;
  vector<string> *gene_names;
  vector<int> *invert_indices;
  vector<long unsigned> *indptr;
  vector<long unsigned> *invert_indptr;
};

void destroy_loaded_data(struct loaded_data_t *loaded_data);
struct loaded_data_t *load_data(const string *path, int use_reduce = 0);
struct loaded_data_t *load_data(const string *path, const string *group_hdf5, int use_reduce = 0);
std::unordered_map<int, int> transform_features_to_go(std::vector<std::string> features,
                              json &gse2indexgo,
                              std::unordered_map<std::string, std::vector<int>> &gse2indexgo_stl,
                              std::unordered_map<std::string, bool> &bad_features);
void write_hdf5(const std::string* path,
                vvvp matrix,
                const std::string* study_name);
void write_matrix_hdf5(std::string group,
                      std::string path,
                      std::vector<float> *data,
                      std::vector<int> *shape,
                      std::vector<long unsigned> *indptr,
                      std::vector<int> *indices,
                      std::vector<std::string> *barcodes,
                      std::vector<std::string> *gene_names);
#endif
