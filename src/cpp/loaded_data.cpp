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

struct loaded_data_t *load_data(const string *path, int use_reduce)
{
	auto res = (struct loaded_data_t *) malloc(sizeof(struct loaded_data_t));
	File file(*path);
	Group group;
	group = file.getGroup("bioturing");
	fprintf(stdout, "Loading indptr\n");
	res->indptr = extract1DArray<long unsigned>(group, "indptr");
	fprintf(stdout, "Loading indices\n");
	res->indices = extract1DArray<int>(group, "indices");
  fprintf(stdout, "Loading data\n");
  res->data = extract1DArray<float>(group, "data");

	if (use_reduce) {
		int n_cells = int(res->invert_indptr->size() - 1);
		int n_genes = int(res->indptr->size() - 1);
		res->shape = new vector<int>({n_cells, n_genes});
		res->gene_names = new vector<string>(n_genes);
		for (int i = 0; i < n_genes; ++i)
			(*(res->gene_names))[i] = to_string(i);
	} else {
		Group group = file.getGroup("countsT");
		res->shape = extract1DArray<int>(group, "shape");
		res->gene_names = extract1DArray(group, "features");
	}

  return res;
}

void write_hdf5(const std::string* path,
                vvvp matrix,
                const std::string* study_name)
{
  File file(pathAppend(*path, "matrix.hdf5"), File::ReadWrite | File::Create | File::Truncate);
  std::vector<std::string> barcodes;
  std::vector<int> indptr = {0};
  std::vector<int> indices;
  std::vector<float> data;
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

  Group group = file.getGroup("/");
  write1DArray(group, "features", list_gene);
  write1DArray(group, "barcodes", barcodes);
  write1DArray(group, "indices", indices);
  write1DArray(group, "indptr", indptr);
  write1DArray(group, "shape", shape);
  write1DArray<float>(group, "data", data);
}

void destroy_loaded_data(struct loaded_data_t *loaded_data)
{
	delete loaded_data->indptr;
	delete loaded_data->indices;
	delete loaded_data->gene_names;
	delete loaded_data->shape;
  delete loaded_data->data;
}


