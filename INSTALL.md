Compiling and installing
========================

Cmake installation
------------------

libNEGF can be installed using `cmake`. The suggested way of installing
is performing an out-of-source build in a temporary directory, e.g. from
the base directory:

```
mkdir build
cd build
cmake ..
make && make install
```

For a list of available configuration options type

```
cmake .. -LAH
```

If you want to build also the tests, enable the `BUILD_TESTING` option:

```
cmake .. -DBUILD_TESTING=TRUE
make && make test
```

MPI support
------------

libNEGF can be compiled with parallel support. Note that in this case
the library has to be linked to [mpifx](https://github.com/dftbplus/mpifx),
which needs to be built separately.

MPI support is enabled with the cmake flag `-DWITH_MPI=TRUE`.