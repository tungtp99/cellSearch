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

float* get_mean_rank(const std::string* path_hdf5, const int* cluster_id, int number_cluster_cells) {
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

  std::cout << " NEW Get mean rank \n";

  float count = 0;
  int number_divide = number_cluster_cells;
  --number_cluster_cells;
  const int* cluster_iter = cluster_id;

  for (int i = 0; i < num_cells; ++i) {
    if (*cluster_iter < i) {
      if (number_cluster_cells == 0) break;
      cluster_iter++;
      --number_cluster_cells;
    }
    if (*cluster_iter != i) {
        pos += (*indptr)[i + 1] - (*indptr)[i];
        continue;
    }
    // printProgress(1.0 * i / (num_cells - 1));

    for (int j = (*indptr)[i]; j < (*indptr)[i + 1]; ++j) {
      row_data.push_back(std::make_pair((*data)[pos], (*indices)[pos]));
      pos++;
    }

    int num_zeros = num_genes - row_data.size();
    if (row_data.size()) {
      std::sort(row_data.begin(), row_data.end());

      rankking[0] = 1;
      for (int j = 1; j < row_data.size(); ++j)
        if (row_data[j].first == row_data[j - 1].first)
          rankking[j] = rankking[j - 1] + 1;
        else 
          rankking[j] = 1;

      for (int j = row_data.size() - 2; j >= 0; --j)
        if (row_data[j].first == row_data[j + 1].first)
          rankking[j] = rankking[j + 1];
      
      count += float(0.5) * (num_genes - row_data.size()) / number_divide;
      int count_smaller = num_genes - row_data.size();
      float rank_current = double(0.5) * (rankking[0] + 1) + count_smaller;
      //float rank_current = 0;
        
      std::map<int, int> check;

      mean_rank[row_data[0].second] += rank_current;
      for (int j = 1; j < row_data.size(); ++j) {
        if (row_data[j].first != row_data[j - 1].first) {
          count_smaller += rankking[j - 1];
          rank_current = double(0.5) * (rankking[j] + 1) + count_smaller;
          if (rank_current > num_genes) {
              std::cout <<"Some thing wrong " << count_smaller << rankking[j] << "\n";
          }
        }
        if (check[row_data[j].second] != 0)  {
            std::cout << "cols appear more then two times" << " " <<  row_data[j].first << " " << row_data[j].second << "\n";
        }
        check[row_data[j].second] = 1;
        mean_rank[row_data[j].second] += rank_current;// - float(0.5) * (num_genes - row_data.size());
      }
      row_data.clear();
    }
  }
  std::cout << "\n";


  std::cout << num_genes;
  for (int i = 0; i < num_genes; i++) {
    mean_rank[i] /= number_divide;
    //mean_rank[i] += count;
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


// Get data of cluster from matrix
vvp filter_cluster_data(loaded_data_t* loaded_data,
                        const int* cluster_ids,
                        const int cluster)
{
  std::cout << "Filter cluster data for cluster " << cluster << "\n";
  vvp cluster_data;
  std::vector<long unsigned int> *indptr = loaded_data->indptr;
  std::vector<int> *indices = loaded_data->indices;
  std::vector<float> *data = loaded_data->data;
  std::vector<int> *shape = loaded_data->shape;
  int num_cells = (*shape)[0];
  int num_genes = (*shape)[1];
  int pos = 0;
  for (int i = 0; i < num_cells; ++i) {
    // std::cout << cluster_ids[i] << " " << cluster << "\n";
    if (cluster_ids[i] != cluster) {
      pos += (*indptr)[i + 1] - (*indptr)[i];
      continue;
    }
    vp cell_data;
    for (int j = (*indptr)[i]; j < (*indptr)[i + 1]; ++j) {
      cell_data.push_back(std::make_pair((*data)[pos],
                                        (*indices)[pos]));
      pos++;
    }
    cluster_data.push_back(std::move(cell_data));
  }
  return std::move(cluster_data);
}

// Represent a cluster by 100 cells which 
// have genes express as mean of 100 random 
// cells in this cluster
vvp make_represent_cluster(loaded_data_t* loaded_data, 
                            const int number_genes,
                            const int* cluster_ids,
                            const int cluster)
{
  std::vector<int> v(100);
  float express_per_row[number_genes];
  std::memset(express_per_row, 0, number_genes * sizeof(float));

  vvp cluster_data = filter_cluster_data(loaded_data,
                                        cluster_ids,
                                        cluster);
  float range =  (1.0 * RAND_MAX + 1u) / cluster_data.size();
  // Get some mean cell represent for cluster 
  vvp cluster_represent;
  std::cout << "cluster " << cluster << " data have size " <<  cluster_data.size() << "\n";

  if (cluster_data.size() == 0)
    return cluster_represent;

  for (int i = 0; i < 100; ++i) {
    
    vp cell_represent;
    // Get random some cell in cluster to calc mean express
    std::generate(v.begin(), v.end(), std::rand);

    for (int j = 0; j < 100; ++j) {
      int cell_index = v[j] / range;
      vp *cell_data = &cluster_data[cell_index];
      for (int k = 0; k < cell_data->size(); k++) 
        express_per_row[(*cell_data)[k].second] += (*cell_data)[k].first;
    }
    // Calc mean express for 100 cells
    for (int j = 0; j < 100; ++j) {
      int cell_index = v[i] / range;
      vp *cell_data = &cluster_data[cell_index];
      for (int k = 0; k < cell_data->size(); k++) {
        int pos = (*cell_data)[k].second;
        if (express_per_row[pos] != 0) {
          cell_represent.push_back(std::make_pair(express_per_row[pos] / 100,
                                                  pos));
          express_per_row[pos] = 0;
        }
      }
    }
    cluster_represent.push_back(std::move(cell_represent));
  }
  
  return std::move(cluster_represent);
}

// For each cluster, calc a group represent for this cluster
// by 100 cells 
void make_represent_matrix(const std::string* path_hdf5,
                          const std::string* path_represent_hdf5,
                          int number_cells,
                          int* cluster_ids,
                          int number_clusters) 
{
  struct loaded_data_t* loaded_data = load_data(path_hdf5);
  std::vector<int> *shape = loaded_data->shape;
  int num_cells = (*shape)[0];
  int num_genes = (*shape)[1];
  if (number_cells != num_cells) {
    std::cout << "Error!!! Number of cells is different with matrix" << endl;
    return;
  }
  vvp cluster_data;
  for (int i = 0; i < number_clusters; ++i) {
    std::cout << "Make represent data for cluster " << i << "of " <<  number_clusters<< "\n";
    vvp represent_cluster = make_represent_cluster(loaded_data,
                                                  num_genes,
                                                  cluster_ids,
                                                  i);
    std::cout << "Represent size" << represent_cluster.size() << "\n";
    cluster_data.insert(cluster_data.end(),
                        represent_cluster.begin(),
                        represent_cluster.end());
  }
  
  write_hdf5(path_represent_hdf5, std::move(cluster_data));
  destroy_loaded_data(loaded_data);
}

extern "C" 
  {
    float* c_get_mean_rank(char* path_hdf5,
                          int* cluster_id,
                          int number_cells) 
    {
      std::string path(path_hdf5);
      return get_mean_rank(&path, cluster_id, number_cells);
    }

    void c_make_represent_matrix(char* path_hdf5_char,
                          char* path_represent_hdf5_char,
                          int number_cells,
                          int* cluster_ids,
                          int number_clusters)
    {
      std::string path_hdf5(path_hdf5_char);
      std::string path_represent_hdf5(path_represent_hdf5_char);
      make_represent_matrix(&path_hdf5,
                            &path_represent_hdf5,
                            number_cells,
                            cluster_ids,
                            number_clusters);
    }

    void c_free_mem(float* ptr) {
      delete[] ptr;
    }
  }