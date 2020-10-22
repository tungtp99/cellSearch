#ifndef __LOADED_DATA__
#define __LOADED_DATA__
#include <vector>
#include <string>
#include "json.hpp"


typedef std::vector<std::pair<float, int>> vp;
typedef std::vector<std::vector<std::pair<float, int>>> vvp;

using namespace std;
using json = nlohmann::json;

struct loaded_data_t {
	vector<long unsigned> *indptr;
	vector<int> *indices;
	vector<long unsigned> *invert_indptr;
	vector<int> *invert_indices;
	vector<string> *gene_names;
	vector<int> *shape;
  vector<float> *data;
	json *annotation;
};

void destroy_loaded_data(struct loaded_data_t *loaded_data);
struct loaded_data_t *load_data(const string *path, int use_reduce = 0);
void write_hdf5(const std::string* path, vvp matrix);
#endif
