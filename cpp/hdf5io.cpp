#include "hdf5io.h"
#include <string.h>

template<typename T>
vector<T> *extract1DArray(const Group &group, const char *name)
{
	DataSet dataset = group.getDataSet(name);
	vector<T> *res = new vector<T>();
	dataset.read(*res);
	return res;
}

vector<string> *extract1DArray(const Group &group, const char *name) {
	DataSet dataset = group.getDataSet(name);
	DataSpace dataSpace = dataset.getSpace();
	auto dims = dataSpace.getDimensions();
	DataType dataType = dataset.getDataType();
	size_t str_size = H5Tget_size(dataType.getId());
	H5T_str_t str_pad = H5Tget_strpad(dataType.getId());
	H5T_cset_t str_cset = H5Tget_cset(dataType.getId());

	vector<string> *res = new vector<string>();
	if((str_pad == H5T_STR_NULLTERM) && (str_cset == H5T_CSET_UTF8)) {
		dataset.read(*res);
	} else {
		dataset.read(*res, str_size, str_pad, str_cset, dims[0]);
	}
	for (long unsigned i = 0; i < res->size(); ++i)
		(*res)[i] = string((*res)[i].c_str());
	return res;
}

template<typename T>
void write1DArray(Group &group, const char *name, const vector<T> &data)
{
	DataSet dataset = group.createDataSet<T>(name, DataSpace::From(data));
	dataset.write(data);
}

template<typename T>
vector<T> *extractRange1DArray(const File &file, const string &groupName,
				const string &datasetName,
				long unsigned start,
				long unsigned end)
{
	DataSet dataset = file.getDataSet(groupName + "/" + datasetName);
	auto res = new vector<T>();
	dataset.select({start}, {end - start}).read(*res);
	return res;
}

vector<string> *extractRange1DArray(const File &file, const string &groupName,
					const std::string &datasetName,
					long unsigned start,
					long unsigned end)
{
	DataSet dataset = file.getDataSet(groupName + "/" + datasetName);
	DataSpace space = dataset.getSpace();
	auto dims = space.getDimensions();
	DataType dataType = dataset.getDataType();
	size_t str_size = H5Tget_size(dataType.getId());
	H5T_str_t str_pad = H5Tget_strpad(dataType.getId());
	H5T_cset_t str_cset = H5Tget_cset(dataType.getId());
	auto res = new vector<string>();
	dataset.select({start}, {end - start}).read(*res, str_size, str_pad,
							str_cset, dims[0]);
	return res;
}

template vector<int> *extract1DArray<int>(const Group &, const char *);
template vector<float> *extract1DArray<float>(const Group &, const char *);
template vector<long unsigned> *extract1DArray<long unsigned>(const Group &, const char *);
template void write1DArray<int>(Group &, const char *, const vector<int> &);
template void write1DArray<float>(Group &, const char *, const vector<float> &);
template void write1DArray<long unsigned>(Group &, const char *, const vector<long unsigned> &);
template void write1DArray<string>(Group &, const char *, const vector<string> &);
template vector<int> *extractRange1DArray(const File &, const string &, const string &,
						long unsigned,
						long unsigned);
template vector<float> *extractRange1DArray(const File &, const string &, const string &,
						long unsigned,
						long unsigned);
template vector<long unsigned> *extractRange1DArray(const File &, const string &, const string &,
						long unsigned,
						long unsigned);

