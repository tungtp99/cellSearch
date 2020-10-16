#include "loaded_data.h"
#include "hdf5io.h"

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

void destroy_loaded_data(struct loaded_data_t *loaded_data)
{
	delete loaded_data->indptr;
	delete loaded_data->indices;
	delete loaded_data->gene_names;
	delete loaded_data->shape;
  delete loaded_data->data;
}


