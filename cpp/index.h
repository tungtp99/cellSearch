#ifndef __INDEX__
#define __INDEX__

#include <iostream>
#include <string>
#include <highfive/H5File.hpp>
#include "loaded_data.h"
#include <vector>
#include <map>

float* get_mean_rank(const std::string* path_hdf5);

#endif