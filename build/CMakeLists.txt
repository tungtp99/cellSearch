tung:
    g++ -shared -Wl,-soname,libindex.so -o libindex.so -fPIC index.cpp hdf5io.cpp loaded_data.cpp -I/usr/local/include -I/home/tung/signac/inst `pkg-config --libs --cflags hdf5`
