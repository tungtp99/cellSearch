#include "index.h"
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

float* get_mean_rank(const std::string* path_hdf5) {
  struct loaded_data_t *loaded_data = load_data(path_hdf5);

  // Get ranking for each row
  std::vector<std::pair<float, int> > row_data;
  std::vector<int> *shape = loaded_data->shape;
  std::vector<long unsigned int> *indptr = loaded_data->indptr;
  std::vector<int> *indices = loaded_data->indices;
  std::vector<float> *data = loaded_data->data;

  int pos = 0;
  int total_add = 0;
  int num_cells = (*shape)[0];
  int num_genes = (*shape)[1];
  float* mean_rank = new float[num_genes];
  std::memset(mean_rank, 0, num_genes * sizeof(float));
  int* rankking = new int[num_genes];

  std::cout << "Get mean rank \n";
  int count = 0;
  for (int i = 0; i < num_cells; ++i) {
    printProgress(1.0 * i / (num_cells - 1));

    for (int j = (*indptr)[i]; j < (*indptr)[i + 1]; ++j) {
      row_data.push_back(std::make_pair((*data)[pos], (*indices)[pos]));
      pos++;
    }

    int num_zeros = num_genes - row_data.size();
    if (row_data.size()) {
      std::sort(row_data.begin(), row_data.end());

      rankking[0] = 0;
      for (int j = 1; j < row_data.size(); ++j)
        if (row_data[j].first == row_data[j - 1].first)
          rankking[j] = rankking[j - 1] + 1;
        else 
          rankking[j] = 1;

      for (int j = row_data.size() - 2; j >= 0; --j)
        if (row_data[j].first == row_data[j + 1].first)
          rankking[j] = rankking[j + 1];
      
      int count_smaller = 0;
      float rank_current = 0.5 * rankking[0] * (rankking[0] + 1);
      mean_rank[row_data[0].second] += rank_current;
      for (int j = 1; j < row_data.size(); ++j) {
        if (row_data[j].first != row_data[j - 1].first) {
          count_smaller += rankking[j - 1];
          rank_current = 0.5 * rankking[j] * (rankking[j] + 1) + count_smaller;
        }
        mean_rank[row_data[j].second] += rank_current;
      }
      row_data.clear();
    }
  }
  std::cout << "\n";

  for (int i = 0; i < num_genes; i++) {
    mean_rank[i] /= num_cells;
  }

  std::cout << "Some element of mean rank ..." << endl;
  for (int i = 0; i < min(10, num_cells); ++i) {
    std::cout << mean_rank[i] << " ";
  }
  std::cout << "\n";


  delete[] rankking;
  destroy_loaded_data(loaded_data);
  return mean_rank;
}


extern "C" {
    float* c_get_mean_rank(char* path_hdf5) {
      std::string path(path_hdf5);
      return get_mean_rank(&path);
    }

    void c_free_mem(float* ptr) {
      delete[] ptr;
    }
  }