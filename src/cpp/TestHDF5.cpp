#include <iostream>
#include <string>
#include <highfive/H5File.hpp>
#include "loaded_data.h"
#include "index.h"

int main (void)
{
  std::string a = "/home/tung/.BioTBDataDev/Data/SingleCell/Study/TabulaMuris_droplet/main/matrix.hdf5";
  int x = 1;
  get_mean_rank(&a);
  return 0;
}
