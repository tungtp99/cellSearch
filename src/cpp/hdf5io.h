#ifndef __HDF5IO__
#define __HDF5IO__
#include <highfive/H5File.hpp>
#include <vector>
#include <string>
using namespace std;
using namespace HighFive;

template<typename T>
//vector<T> *extract1DArray(const Group &group, const char *name,
//				const PredType &predtype);
vector<T> *extract1DArray(const Group &group, const char *name);

vector<string> *extract1DArray(const Group &group, const char *name);

template<typename T>
//void write1DArray(Group &group, const char *name, const vector<T> &data,
//			const PredType &predtype);
void write1DArray(Group &group, const char *name, const vector<T> &data);

template<typename T>
vector<T> *extractRange1DArray(const File &file, const string &groupName,
				const string &datasetName,
				long unsigned start,
				long unsigned end);

vector<string> *extractRange1DArray(const File &file, const string &groupName,
					const std::string &datasetName,
					long unsigned start,
					long unsigned end);
#endif
